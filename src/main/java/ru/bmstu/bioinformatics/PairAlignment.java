package ru.bmstu.bioinformatics;

import java.util.Arrays;

public class PairAlignment {
    private final String OUTPUT_SEQUENCE_1 = "Seq1: ", OUTPUT_SEQUENCE_2 = "Seq2: ", OUTPUT_SCORE = "Score: ",
            NEXT_LINE = "\n";
    private final char GAP = '_';
    private final int ZERO = 0, GAP_INDEX = 0, ONE = 1, FIRST_INDEX = 1, NUMBER_OF_SYMBOLS_IN_LINE = 50;

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

        for (int columnIndex = FIRST_INDEX; columnIndex <= columnSize; columnIndex++) {
            currentMaximum = matrix[lineSize][columnIndex].getValue();
            if (currentMaximum > maximum) {
                columnMaximum = columnIndex;
                maximum = currentMaximum;
            }
        }
        for (int lineIndex = FIRST_INDEX; lineIndex <= lineSize; lineIndex++) {
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
        char firstChar = firstSequence.charAt(columnIndex - ONE),
                secondChar = secondSequence.charAt(lineIndex - ONE);
        if (firstChar == secondChar) {
            return getMatch(firstChar, secondChar);
        } else return getMismatch(firstChar, secondChar);
    }

    private void fillGapCells() {
        matrix[GAP_INDEX][GAP_INDEX] = new Cell(null, PredecessorType.NULL, GAP, GAP, ZERO);

        for (int columnIndex = FIRST_INDEX; columnIndex < matrix[GAP_INDEX].length; columnIndex++) {
            int value = optimization
                    ? ZERO
                    : getIndel() * columnIndex;
            matrix[GAP_INDEX][columnIndex] = new Cell(
                    matrix[GAP_INDEX][columnIndex - ONE],
                    PredecessorType.LEFT,
                    firstSequence.charAt(columnIndex - ONE),
                    GAP,
                    value);
        }
        for (int lineIndex = FIRST_INDEX; lineIndex < matrix.length; lineIndex++) {
            int value = optimization
                    ? ZERO
                    : getIndel() * lineIndex;
            matrix[lineIndex][GAP_INDEX] = new Cell(
                    matrix[lineIndex - ONE][GAP_INDEX],
                    PredecessorType.UP,
                    GAP,
                    secondSequence.charAt(lineIndex - ONE),
                    value);
        }
    }

    private void fillScoringMatrix() {
        for (int lineIndex = FIRST_INDEX; lineIndex < matrix.length; lineIndex++) {
            for (int columnIndex = FIRST_INDEX; columnIndex < matrix[GAP_INDEX].length; columnIndex++) {
                int currentMaximum = matrix[lineIndex - ONE][columnIndex - ONE].getValue()
                        + checkMatch(lineIndex, columnIndex);
                PredecessorType currentType = PredecessorType.DIAG;

                if (matrix[lineIndex - 1][columnIndex].getValue() + getIndel() > currentMaximum) {
                    currentMaximum = matrix[lineIndex - ONE][columnIndex].getValue() + getIndel();
                    currentType = PredecessorType.UP;
                }
                if (matrix[lineIndex][columnIndex - ONE].getValue() + getIndel() > currentMaximum) {
                    currentMaximum = matrix[lineIndex][columnIndex - ONE].getValue() + getIndel();
                    currentType = PredecessorType.LEFT;
                }

                Cell predecessorCell = null;
                switch (currentType) {
                    case LEFT:
                        predecessorCell = matrix[lineIndex][columnIndex - ONE];
                        break;
                    case DIAG:
                        predecessorCell = matrix[lineIndex - ONE][columnIndex - ONE];
                        break;
                    case UP:
                        predecessorCell = matrix[lineIndex - ONE][columnIndex];
                        break;
                }

                matrix[lineIndex][columnIndex] = new Cell(
                        predecessorCell,
                        currentType,
                        firstSequence.charAt(columnIndex - ONE),
                        secondSequence.charAt(lineIndex - ONE),
                        currentMaximum);
            }
        }
    }

    private void align() {
        matrix = new Cell[secondSequence.length() + ONE][firstSequence.length() + ONE];

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
        int i = FIRST_INDEX;
        for (; i < Math.ceil((double) firstAlignSequence.length() / NUMBER_OF_SYMBOLS_IN_LINE); i++) {
            stringBuilder.append(OUTPUT_SEQUENCE_1);
            stringBuilder.append(firstAlignSequence.substring(
                    NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE),
                    NUMBER_OF_SYMBOLS_IN_LINE * i));
            stringBuilder.append(NEXT_LINE);

            stringBuilder.append(OUTPUT_SEQUENCE_2);
            stringBuilder.append(secondAlignSequence.substring(
                    NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE),
                    NUMBER_OF_SYMBOLS_IN_LINE * i));
            stringBuilder.append(NEXT_LINE);
            stringBuilder.append(NEXT_LINE);
        }

        stringBuilder.append(OUTPUT_SEQUENCE_1);
        stringBuilder.append(firstAlignSequence.substring(NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE)));
        stringBuilder.append(NEXT_LINE);

        stringBuilder.append(OUTPUT_SEQUENCE_2);
        stringBuilder.append(secondAlignSequence.substring(NUMBER_OF_SYMBOLS_IN_LINE * (i - ONE)));
        stringBuilder.append(NEXT_LINE);
        stringBuilder.append(NEXT_LINE);

        stringBuilder.append(OUTPUT_SCORE);
        stringBuilder.append(score);

        return stringBuilder.toString();
    }
}
