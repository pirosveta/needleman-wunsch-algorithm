package ru.bmstu.bioinformatics;

import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

@RunWith(Parameterized.class)
public class PairAlignmentTest {

    @Parameterized.Parameters(name = "Test {index}")
    public static Iterable<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {"./seq1_1.fasta", "./seq2_1.fasta", "BLOSUM62", "-1", 293, 312},
                {"./seq1_2.fasta", "./seq2_2.fasta", "DNAFull", "-5", -5, 15},
                {"./seq1_3.fasta", "./seq2_3.fasta", "Default", "-5", -19, 5}
        });
    }

    private final String SPACE_REGEX = "\\s+", OUTPUT_FILE_PATH = "./out.txt", EMPTY = "",
        PATH_KEY = "-i", ALPHABET_KEY = "-a", GAP_PENALTY_KEY = "-g", OUTPUT_FILE_KEY = "-o",
            OPTIMIZATION_KEY = "-optimization", TRUE = "true", FALSE = "false";
    private final int ZERO = 0, SCORE_INDEX = 1, EXIT_CODE_WITH_ERROR = 1;

    private String firstSequencePath,
            secondSequencePath,
            alphabet,
            gapPenalty;
    private int notOptimizedScore, optimizedScore;

    private int getScore() {
        int score = ZERO;

        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(OUTPUT_FILE_PATH));
            String currentLine = EMPTY, nextLine;

            while ((nextLine = bufferedReader.readLine()) != null) {
                currentLine = nextLine;
            }
            String[] splitLine = currentLine.split(SPACE_REGEX);
            score = Integer.parseInt(splitLine[SCORE_INDEX]);

            bufferedReader.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(EXIT_CODE_WITH_ERROR);
        }

        return score;
    }

    public PairAlignmentTest(String firstSequencePath, String secondSequencePath,
                             String alphabet, String gapPenalty, int notOptimizedScore, int optimizedScore) {
        this.firstSequencePath = firstSequencePath;
        this.secondSequencePath = secondSequencePath;
        this.alphabet = alphabet;
        this.gapPenalty = gapPenalty;
        this.notOptimizedScore = notOptimizedScore;
        this.optimizedScore = optimizedScore;
    }

    @Test
    public void test1() {
        Main.main(new String[]{
                PATH_KEY, firstSequencePath, secondSequencePath,
                ALPHABET_KEY, alphabet,
                GAP_PENALTY_KEY, gapPenalty,
                OUTPUT_FILE_KEY, OUTPUT_FILE_PATH,
                OPTIMIZATION_KEY, FALSE});
        Assert.assertEquals(getScore(), notOptimizedScore);
    }

    @Test
    public void test2() {
        Main.main(new String[]{
                PATH_KEY, firstSequencePath, secondSequencePath,
                ALPHABET_KEY, alphabet,
                GAP_PENALTY_KEY, gapPenalty,
                OUTPUT_FILE_KEY, OUTPUT_FILE_PATH,
                OPTIMIZATION_KEY, TRUE});
        Assert.assertEquals(getScore(), optimizedScore);
    }
}
