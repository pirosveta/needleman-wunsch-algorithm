package ru.bmstu.bioinformatics;

public class Cell {
    private final char GAP = '_';
    private final int BEGIN = 0;

    private Cell predecessor;
    private PredecessorType predecessorType;

    private char inFirstSequence;
    private char inSecondSequence;
    private int value;

    public Cell(Cell predecessor, PredecessorType predecessorType,
                char inFirstSequence, char inSecondSequence,
                int value) {
        this.predecessor = predecessor;
        this.predecessorType = predecessorType;
        this.inFirstSequence = inFirstSequence;
        this.inSecondSequence = inSecondSequence;
        this.value = value;
    }

    public Cell fillInformation(StringBuilder firstSequence, StringBuilder secondSequence) {
        if (predecessor != null) {
            switch (predecessorType) {
                case LEFT:
                    firstSequence.insert(BEGIN, inFirstSequence);
                    secondSequence.insert(BEGIN, GAP);
                    break;
                case DIAG:
                    firstSequence.insert(BEGIN, inFirstSequence);
                    secondSequence.insert(BEGIN, inSecondSequence);
                    break;
                case UP:
                    firstSequence.insert(BEGIN, GAP);
                    secondSequence.insert(BEGIN, inSecondSequence);
                    break;
            }
        }

        return predecessor;
    }

    public int getValue() {
        return value;
    }
}
