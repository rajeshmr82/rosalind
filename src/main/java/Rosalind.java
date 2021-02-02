package main.java;

import java.io.*;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Objects;
import java.util.regex.Pattern;

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

    // 5. GC 	Computing GC Content
    public static void maxGCContent(String filename) throws IOException {
        File inputFile=new File(Objects.requireNonNull(Solution.class.getClassLoader().getResource(filename)).getFile());
        InputStream inputStream = new FileInputStream(inputFile);
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        String line;
        String maxID = null;
        double maxValue=0.0;
        String id = null;
        HashMap<String,String> map=new HashMap<>();
        while ((line = reader.readLine()) != null) {
            var matcher = Pattern.compile("(>Rosalind_\\d{4})(.*)").matcher(line);
            if(matcher.find()){
                id= matcher.group(1);
                map.put(id,matcher.group(2));
            }else{
                map.put(id,map.getOrDefault(id,"").concat(line));
            }
        }

        for (var v:
             map.entrySet()) {
            var countCG = v.getValue().chars().filter(ch -> ch=='C' || ch=='G').count();
            var gcContent = countCG*100.0/v.getValue().length();
            System.out.println("Id:"+id+" DNA:"+v.getValue()+ " CG:"+countCG+ " Len:"+ v.getValue().length() + " GCContent:"+gcContent);
            if(gcContent>maxValue){
                maxValue = gcContent;
                maxID = v.getKey();
            }
        }

        System.out.println(maxID);
        DecimalFormat df = new DecimalFormat("##.######");
        System.out.println(df.format(maxValue));
    }

    // 6. HAMM 	Counting Point Mutations
    public static int computeHammingDistance(String s, String t) {
        if(s.length()!=t.length()){
            return -1;
        }
        int distance=0;
        for (int i = 0; i < s.length(); i++) {
            if(s.charAt(i)!=t.charAt(i)){
                distance++;
            }
        }

        return distance;
    }

    // 7. IPRB 	Mendel's First Law
    public static String probabilityOfDominant(int k, int m, int n) {
        var total = k + m + n;
        var pM = m * 1.0 / total;
        var pN = n * 1.0 / total;
        DecimalFormat df = new DecimalFormat("#.#####");
        return df.format(1 - pN * ((n - 1.0) / (total - 1.0)) - pN * (m / (total - 1.0)) - pM * ((m - 1.0) / (total - 1.0)) * 0.25);
    }

}