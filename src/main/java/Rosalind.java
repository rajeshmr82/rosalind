package main.java;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

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

    // 8. PROT 	Translating RNA into Protein
    public static String translateRNA(String mRNA) {
        HashMap<String, String> rnaCodon = initRNACodon();

        StringBuilder result = new StringBuilder();
        for (int i = 0; i < mRNA.length(); i=i+3) {
            if (rnaCodon.getOrDefault(mRNA.substring(i, i + 3),"").equals("Stop")) {
                break;
            }
            result.append(rnaCodon.getOrDefault(mRNA.substring(i, i + 3), ""));
        }

        return result.toString();
    }

    private static HashMap<String, String> initRNACodon() {
        HashMap<String, String> rnaCodon = new HashMap<>();
        rnaCodon.put("UUU","F");
        rnaCodon.put("CUU","L");
        rnaCodon.put("AUU","I");
        rnaCodon.put("GUU","V");
        rnaCodon.put("UUC","F");
        rnaCodon.put("CUC","L");
        rnaCodon.put("AUC","I");
        rnaCodon.put("GUC","V");
        rnaCodon.put("UUA","L");
        rnaCodon.put("CUA","L");
        rnaCodon.put("AUA","I");
        rnaCodon.put("GUA","V");
        rnaCodon.put("UUG","L");
        rnaCodon.put("CUG","L");
        rnaCodon.put("AUG","M");
        rnaCodon.put("GUG","V");
        rnaCodon.put("UCU","S");
        rnaCodon.put("CCU","P");
        rnaCodon.put("ACU","T");
        rnaCodon.put("GCU","A");
        rnaCodon.put("UCC","S");
        rnaCodon.put("CCC","P");
        rnaCodon.put("ACC","T");
        rnaCodon.put("GCC","A");
        rnaCodon.put("UCA","S");
        rnaCodon.put("CCA","P");
        rnaCodon.put("ACA","T");
        rnaCodon.put("GCA","A");
        rnaCodon.put("UCG","S");
        rnaCodon.put("CCG","P");
        rnaCodon.put("ACG","T");
        rnaCodon.put("GCG","A");
        rnaCodon.put("UAU","Y");
        rnaCodon.put("CAU","H");
        rnaCodon.put("AAU","N");
        rnaCodon.put("GAU","D");
        rnaCodon.put("UAC","Y");
        rnaCodon.put("CAC","H");
        rnaCodon.put("AAC","N");
        rnaCodon.put("GAC","D");
        rnaCodon.put("CAA","Q");
        rnaCodon.put("AAA","K");
        rnaCodon.put("GAA","E");
        rnaCodon.put("CAG","Q");
        rnaCodon.put("AAG","K");
        rnaCodon.put("GAG","E");
        rnaCodon.put("UGU","C");
        rnaCodon.put("CGU","R");
        rnaCodon.put("AGU","S");
        rnaCodon.put("GGU","G");
        rnaCodon.put("UGC","C");
        rnaCodon.put("CGC","R");
        rnaCodon.put("AGC","S");
        rnaCodon.put("GGC","G");
        rnaCodon.put("CGA","R");
        rnaCodon.put("AGA","R");
        rnaCodon.put("GGA","G");
        rnaCodon.put("UGG","W");
        rnaCodon.put("CGG","R");
        rnaCodon.put("AGG","R");
        rnaCodon.put("GGG","G");
        rnaCodon.put("UAA","Stop");
        rnaCodon.put("UAG","Stop");
        rnaCodon.put("UGA","Stop");
        return rnaCodon;
    }

    // 9. SUBS 	Finding a Motif in DNA
    public static String findAMotif(String s, String t) {
        if(t.length()>s.length()){
            return "0";
        }

        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < s.length()-t.length(); i++) {
            if(s.startsWith(t, i)){
                indices.add(i+1);
            }
        }

        return indices
                .stream()
                .map(String::valueOf)
                .collect(Collectors.joining(" "));
    }
}