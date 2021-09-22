package ru.bmstu.bioinformatics;

import java.util.Arrays;

public class PairAlignment {
    private final char GAP = '_';

    private Cell[][] matrix;
    private String firstSequence,
            secondSequence;
    private ScoringFunction scoringFunction;
    private boolean optimization;

    private StringBuilder firstAlignSequence = new StringBuilder(),
            secondAlignSequence = new StringBuilder();
    private int score;

    private int getIndel() {
        if (scoringFunction instanceof Default) {
            return ((Default) scoringFunction).getIndel();
        } else if (scoringFunction instanceof DNAFull) {
            return ((DNAFull) scoringFunction).getIndel();
        } else return ((BLOSUM62) scoringFunction).getIndel();
    }

    private int getMatch(char firstChar, char secondChar) {
        if (scoringFunction instanceof Default) {
            return ((Default) scoringFunction).getMatch();
        } else if (scoringFunction instanceof DNAFull) {
            return ((DNAFull) scoringFunction).getMatch();
        } else return ((BLOSUM62) scoringFunction).getScore(firstChar, secondChar);
    }

    private int getMismatch(char firstChar, char secondChar) {
        if (scoringFunction instanceof Default) {
            return ((Default) scoringFunction).getMismatch();
        } else if (scoringFunction instanceof DNAFull) {
            return ((DNAFull) scoringFunction).getMismatch();
        } else return ((BLOSUM62) scoringFunction).getScore(firstChar, secondChar);
    }

    private void fillMissingCharacters(int index, boolean inFirstSequence) {
        if (inFirstSequence) {
            firstAlignSequence.append(firstSequence.substring(index));
            char[] addition = new char[firstSequence.length() - index];
            Arrays.fill(addition, GAP);
            secondAlignSequence.append(addition);
        } else {
            char[] addition = new char[secondSequence.length() - index];
            Arrays.fill(addition, GAP);
            firstAlignSequence.append(addition);
            secondAlignSequence.append(secondSequence.substring(index));
        }
    }

    private Cell findOptimalScoreCell() {
        int lineSize = secondSequence.length(),
                columnSize = firstSequence.length(),
                lineMaximum = lineSize,
                columnMaximum = columnSize,
                maximum = matrix[lineMaximum][columnMaximum].getValue(),
                currentMaximum;

        for (int columnIndex = 1; columnIndex <= columnSize; columnIndex++) {
            currentMaximum = matrix[lineSize][columnIndex].getValue();
            if (currentMaximum > maximum) {
                columnMaximum = columnIndex;
                maximum = currentMaximum;
            }
        }
        for (int lineIndex = 1; lineIndex <= lineSize; lineIndex++) {
            currentMaximum = matrix[lineIndex][columnSize].getValue();
            if (currentMaximum > maximum) {
                lineMaximum = lineIndex;
                columnMaximum = columnSize;
                maximum = currentMaximum;
            }
        }

        if (lineMaximum != lineSize) {
            fillMissingCharacters(lineMaximum, false);
        } else if (columnMaximum != columnSize) {
            fillMissingCharacters(columnMaximum, true);
        }

        return matrix[lineMaximum][columnMaximum];
    }

    private int checkMatch(int lineIndex, int columnIndex) {
        char firstChar = firstSequence.charAt(columnIndex - 1),
                secondChar = secondSequence.charAt(lineIndex - 1);
        if (firstChar == secondChar) {
            return getMatch(firstChar, secondChar);
        } else return getMismatch(firstChar, secondChar);
    }

    private void fillGapCells() {
        matrix[0][0] = new Cell(null, PredecessorType.NULL, GAP, GAP, 0);

        for (int columnIndex = 1; columnIndex < matrix[0].length; columnIndex++) {
            int value = optimization
                    ? 0
                    : getIndel() * columnIndex;
            matrix[0][columnIndex] = new Cell(
                    matrix[0][columnIndex - 1],
                    PredecessorType.LEFT,
                    firstSequence.charAt(columnIndex - 1),
                    GAP,
                    value);
        }
        for (int lineIndex = 1; lineIndex < matrix.length; lineIndex++) {
            int value = optimization
                    ? 0
                    : getIndel() * lineIndex;
            matrix[lineIndex][0] = new Cell(
                    matrix[lineIndex - 1][0],
                    PredecessorType.UP,
                    GAP,
                    secondSequence.charAt(lineIndex - 1),
                    value);
        }
    }

    private void fillScoringMatrix() {
        for (int lineIndex = 1; lineIndex < matrix.length; lineIndex++) {
            for (int columnIndex = 1; columnIndex < matrix[0].length; columnIndex++) {
                int currentMaximum = matrix[lineIndex - 1][columnIndex - 1].getValue()
                        + checkMatch(lineIndex, columnIndex);
                PredecessorType currentType = PredecessorType.DIAG;

                if (matrix[lineIndex - 1][columnIndex].getValue() + getIndel() > currentMaximum) {
                    currentMaximum = matrix[lineIndex - 1][columnIndex].getValue() + getIndel();
                    currentType = PredecessorType.UP;
                }
                if (matrix[lineIndex][columnIndex - 1].getValue() + getIndel() > currentMaximum) {
                    currentMaximum = matrix[lineIndex][columnIndex - 1].getValue() + getIndel();
                    currentType = PredecessorType.LEFT;
                }

                Cell predecessorCell = null;
                switch (currentType) {
                    case LEFT:
                        predecessorCell = matrix[lineIndex][columnIndex - 1];
                        break;
                    case DIAG:
                        predecessorCell = matrix[lineIndex - 1][columnIndex - 1];
                        break;
                    case UP:
                        predecessorCell = matrix[lineIndex - 1][columnIndex];
                        break;
                }

                matrix[lineIndex][columnIndex] = new Cell(
                        predecessorCell,
                        currentType,
                        firstSequence.charAt(columnIndex - 1),
                        secondSequence.charAt(lineIndex - 1),
                        currentMaximum);
            }
        }
    }

    private void align() {
        matrix = new Cell[secondSequence.length() + 1][firstSequence.length() + 1];

        fillGapCells();
        fillScoringMatrix();

        Cell optimalScoreCell = optimization
                ? findOptimalScoreCell()
                : matrix[secondSequence.length()][firstSequence.length()];
        score = optimalScoreCell.getValue();
        while (optimalScoreCell != null) {
            optimalScoreCell = optimalScoreCell.fillInformation(firstAlignSequence, secondAlignSequence);
        }
    }

    public PairAlignment(String firstSequence, String secondSequence,
                         ScoringFunction scoringFunction, boolean optimization) {
        this.firstSequence = firstSequence;
        this.secondSequence = secondSequence;
        this.scoringFunction = scoringFunction;
        this.optimization = optimization;

        align();
    }

    @Override
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        int i = 1;
        for (; i < Math.ceil(firstAlignSequence.length() / 50.0); i++) {
            stringBuilder.append("Seq1: " + firstAlignSequence.substring(50 * (i - 1), 50 * i) + "\n");
            stringBuilder.append("Seq2: " + secondAlignSequence.substring(50 * (i - 1), 50 * i) + "\n");
            stringBuilder.append("\n");
        }
        stringBuilder.append("Seq1: " + firstAlignSequence.substring(50 * (i - 1)) + "\n");
        stringBuilder.append("Seq2: " + secondAlignSequence.substring(50 * (i - 1)) + "\n");
        stringBuilder.append("\n");

        stringBuilder.append("Score: " + score);
        return stringBuilder.toString();
    }
}
