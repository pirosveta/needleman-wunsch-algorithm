package ru.bmstu.bioinformatics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Stream;

public class BLOSUM62 implements ScoringFunction {
    private int indel;

    private HashMap<Character, Integer> orderInMatrix = new HashMap<>();
    private ArrayList<ArrayList<Integer>> matrix = new ArrayList<>();

    private void initializeMatrix() {
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader("BLOSUM62.dat"));

            String line = bufferedReader.readLine();
            String[] alphabet = line.split("\\s+");
            int delta = 0;
            if (alphabet[0].equals(""))
                delta = 1;

            for (int order = 0 + delta; order < alphabet.length; order++) {
                orderInMatrix.put(alphabet[order].charAt(0), order - delta);
            }
            while ((line = bufferedReader.readLine()) != null) {
                String[] splitLine = line.split("\\s+");
                ArrayList<Integer> scores = new ArrayList<>();
                for (int index = 1; index < splitLine.length; index++) {
                    scores.add(Integer.parseInt(splitLine[index]));
                }
                matrix.add(scores);
            }

            bufferedReader.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public BLOSUM62(int indel) {
        this.indel = indel;
        initializeMatrix();
    }

    public int getIndel() {
        return indel;
    }

    public int getScore(char firstChar, char secondChar) {
        int firstOrder = orderInMatrix.get(firstChar),
                secondOrder = orderInMatrix.get(secondChar);
        return matrix.get(firstOrder).get(secondOrder);
    }
}
