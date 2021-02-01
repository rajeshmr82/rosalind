package main.java;

import java.util.HashMap;

public class Rosalind {
    // 1. DNA 	Counting DNA Nucleotides
    public static String countDNANucleotides(String input) {
        HashMap<Character, Integer> map = new HashMap<>();
        for (Character c :
                input.toCharArray()) {
            map.put(c, map.getOrDefault(c, 0) + 1);
        }

        return map.get('A') + " " + map.get('C') + " " + map.get('G') + " " + map.get('T');
    }

    // 2. RNA 	Transcribing DNA into RNA
    public static String DNAToRNA(String dna) {
        return dna.replaceAll("T", "U");
    }

    // 3. REVC  Complementing a Strand of DNA
    public static String complementDNA(String dna) {
        StringBuilder result = new StringBuilder();
        HashMap<String, String> map = new HashMap<>();
        map.put("A", "T");
        map.put("T", "A");
        map.put("C", "G");
        map.put("G", "C");

        for (int i = dna.length() - 1; i >= 0; i--) {
            result.append(map.get(dna.substring(i, i + 1)));
        }

        return result.toString();
    }

    // 4. FIB 	Rabbits and Recurrence Relations
    public static long rabbitRecurrence(int n, int k) {
        long[] fib = new long[n + 1];
        fib[0] = 0;
        fib[1] = 1;
        for (int i = 2; i < n + 1; i++) {
            fib[i] = fib[i - 1] + k * fib[i - 2];
        }

        return fib[n];
    }
}