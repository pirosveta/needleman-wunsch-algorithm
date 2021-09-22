package ru.bmstu.bioinformatics;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;

@RunWith(Parameterized.class)
public class PairAlignmentTest {

    @Parameterized.Parameters(name = "Test {index}")
    public static Iterable<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {"./seq1_1.fasta", "./seq2_1.fasta", "BLOSUM62", "-1", "-optimization"},
                {"./seq1_2.fasta", "./seq2_2.fasta", "DNAFull", "-5", "-optimization"}
        });
    }

    private String firstSequencePath,
            secondSequencePath,
            alphabet,
            gapPenalty,
            optimization;

    public PairAlignmentTest(String firstSequencePath, String secondSequencePath,
                             String alphabet, String gapPenalty, String optimization) {
        this.firstSequencePath = firstSequencePath;
        this.secondSequencePath = secondSequencePath;
        this.alphabet = alphabet;
        this.gapPenalty = gapPenalty;
        this.optimization = optimization;
    }

    @Test
    public void test1() {
        Main.main(new String[]{"-i", firstSequencePath, secondSequencePath,
                "-a", alphabet,
                "-g", gapPenalty});
    }

    @Test
    public void test2() {
        Main.main(new String[]{"-i", firstSequencePath, secondSequencePath,
                "-a", alphabet,
                "-g", gapPenalty,
                optimization});
    }
}
