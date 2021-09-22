package ru.bmstu.bioinformatics;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class Main {
    @Parameter(names = "-i", description = "Paths to input sequences", arity = 2, required = true)
    private static List<String> inputPaths = new ArrayList<>();

    @Parameter(names = "-a", description = "Alphabet")
    private static String alphabet = "Default";

    @Parameter(names = {"-g", "--gap"}, description = "Penalty for the gap")
    private static String gapPenalty = "";

    @Parameter(names = "-o", description = "Path to output file")
    private static String outputPath = "";

    @Parameter(names = "-optimization", description = "Enable optimization")
    private static boolean optimization = false;

    public static ArrayList<String> readFile() {
        ArrayList<String> sequences = new ArrayList<>();

        for (String sequencePath : inputPaths) {
            try {
                BufferedReader bufferedReader = new BufferedReader(new FileReader(sequencePath));
                StringBuilder stringBuilder = new StringBuilder();
                String line;

                bufferedReader.readLine();
                while ((line = bufferedReader.readLine()) != null) {
                    stringBuilder.append(line);
                }
                sequences.add(stringBuilder.toString());

                bufferedReader.close();
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        return sequences;
    }

    private static ScoringFunction defineScoringFunction() {
        ScoringFunction scoringFunction;

        if (alphabet.equals("DNAFull")) {
            scoringFunction = new DNAFull(Integer.parseInt(gapPenalty), 5, -4);
        }
        else if (alphabet.equals("BLOSUM62")) {
            scoringFunction = new BLOSUM62(Integer.parseInt(gapPenalty));
        }
        else scoringFunction = new Default(2, 1, -1);

        return scoringFunction;
    }

    public static void main(String[] args) {
        try {
            Main main = new Main();
            JCommander jCommander = new JCommander(main);
            jCommander.parse(args);

            if (inputPaths.size() != 2
                    || !(alphabet.equals("BLOSUM62") || alphabet.equals("DNAFull") || alphabet.equals("Default"))
                    || (!alphabet.equals("Default") && gapPenalty.equals(""))) {
                jCommander.usage();
                return;
            }

            ArrayList<String> sequences = readFile();
            PairAlignment pairAlignment = new PairAlignment(
                    sequences.get(0),
                    sequences.get(1),
                    defineScoringFunction(),
                    optimization);

            PrintStream console = System.out;
            if (!outputPath.equals("")) {
                File file = new File(outputPath);
                FileOutputStream fos = new FileOutputStream(file);
                PrintStream ps = new PrintStream(fos);
                System.setOut(ps);
            }
            System.out.println(pairAlignment.toString());
            System.setOut(console);

        } catch (ParameterException | FileNotFoundException e) {
            e.printStackTrace();
        }

    }
}