package main.java;

import java.io.IOException;

public class Solution {
    public static void main(String[] args) throws IOException {
        // 1. DNA 	Counting DNA Nucleotides
        // System.out.println(Rosalind.countDNANucleotides("AACTGAAACATTGCAGGGTTTCACGGGATCATTAATGGCTTCCCGAACAACGCCATCCCCACCAAGTGTGCCAGCCAATGTACAATCCCTGTCATCTAGACGAAACGAGCTATTGTTCATACTGACTTAGTAGTTTCGTACTCCAGTGCAGTAGCCCGCGGGGGGGCGTAGCGTGTGCGGCCGTCAATGGACAATGGATGGCTCGTAGGTGCCAGCAAGCGGAGTCATTTTATACCCCAGCAATTCGACATAGCCTGACACGTTACCGAACCAACCGTCTTCGAGTTTTGAGGAGATAATGAAGTATTAGTGGGGCCAAGGTGGCTTATACGTGTCAAGATCTTGGGGTCCCCAATTAAAGCTCCGTCAACTGACTCATGCACTACTAAATCGATTGGGGATCGTCGTCTGTGTGGTGACTACGAACGGCGGCTCATTTCCAGCGACACGTGCGAAGTACAGGCTGATTCCGCCATGCTCCCGGGTAGCAATTCGTAGATCGGCTTTATCTTTATCACGATTCGTCCTATGATGGTCACATTGCGTCCCACAGCTGTATTCTTATTGCTACATGATGTAAAAACGACTACTGGAATCTAGAGTACTCTACTCACTCCTGCGATCAATTTGAAAGGAGATCACCTACGCAACGGTAGGATGCCCGACAGATGTTGCATTAGTACCATTATAATGAGTCACAAGTCCGCTCGTGAGCCTGCTACTGCGAGTGCCTAGCCTGTATGGGCAAGGCCTATTTGTCTTCAGGCTCTGAGGGTTTACTGGGAATCAGAGCCAATCCGCGTTTTGAGCGCTCGACTCCAGTGAGGCATCATGTGGTAAGATCTAGCGACGGTTCACACTCGCTTGCCCATAGCAGCTGTACCGTACCAACCCAAATACCACTCCGCAACCACGCATGCGCAG"));
        // 2. RNA 	Transcribing DNA into RNA
        // System.out.println(Rosalind.DNAToRNA("ATGCAAAAGGGGAGAAGGTATGCTAACTTGACGCACAGTTAAAAATACCGGGCGGCAGGTCACTGGCATTTGATGAGTGTAAAACTGCGGATAGCGAGAGGCAGTCACCGTCCACAAATTACTCTTTTTGCCCCACGGAAAGAACTCGGTAGACATATTCCACCTAACCATCATGACTCTTCTCACTTAATTCCGTGCTACAAGGGACTTACTCTATGCAGTTTTAAAGCTTTTTGGCCAAGCGTGTCGTTTTTTACCCGAGCTTACTGGCTAGTGCCCCAGTTGGGGTGTAGATGGTATCAGATTAGGTCCGTGAGACTGACTAGCTGCGGGGGCCGACATTGACACGACCCTGGTAGTTGGGGGGCCTAAGACGGGTAGTTATGGTACTTAATCCGCTGCCGGCACATGCGTGAAGTGGCCTCACACAATGCCAACTTCAGTGGTGCCTTCTGCGTCGCACACGACTGTTCGGGACTATCGTGTTAGACAGGATGTACTGAGTTTCAAATACACGTACCTTGTTTCTTGAGTTGGTGCCAGGAGTTGAACATATCAATTGCCCAATCCGCAGTGGGATGTTAGCGCTACTCCTCGTGCTGGCGAAAGAGGCTCGGTATCCTCTGCAGACCTGTTTCAAGGATCGGACAGACAGCGTGGCTTTACCATGATTTCTGTTATTGTTCGACGTTAAGGAAGTATGACGCGCCCAAGTGGACACCTACACGTGTTTCCAAACCTCTACGTCGACGGGAGTAGACAGGGGTGTGCCACAGTATCTTGCCAATAAAAACGCAATCGTCCCTAATTAGCTTTACGAATCGGACATTGTGTGTTTCCGTATCCGACGATATATCGGCATCCCGGGTAGCTACGGTGTCTTTTTCATGTGAGACTGCGTAACACTATTGTTCGCAACAGGGTGCGCACCGCGTGTGCGAAACTTA"));
        // 3. REVC  Complementing a Strand of DNA
        // System.out.println(Rosalind.complementDNA("CGGGGGCGTTGTTCAAACTTTACACTCGACTTTTTGTGGTCACAACTTCGACACGAGTGTACACGGATCAGGGACTTGGTATAATGGGCGGACATCAAAGTCGGCGTATCATGTTCCGCGTGCCAAACACCGGATGCCCATCACAACCGCACGGGATCGTGCGCGACTGCTGAGTGCTTGGTTGCTGCACATTGATTTGATGACTGGGGTTCAGGTTGACTTTTATCGCTTTACGTAGCAAGATATTTTATCACATATGTGCCGTAACCAGTACACGGGACCGGATGACCAGAGACGCGAGTTGGCGACCAAATAGGCCTCATGACAAGAGTAGCCTCGAGACAGCATCTCCCAGTACCTAACCGCTAACCATGTTCTCTTTGCACAATGTAACTGTCGACTCCCGTTACGGGATAATCCCAGCCGCCCGCGATCGGCTAGTCAGTGCCGCCCGTGGAGTGGCACGTGTTACTGAAAGACTCGGTATCATAAATAAGAACCGTGGTGGTCTATAGATGCTAACCAGAGGAAATGTGCTTCTCTTCTTCGTCCCTGGCCCATCTAAACTAGTCGTGGTATCACAAGCATAATCTTCGCGTTCTGGACTGTACCCAAGCCACCTACCGCGACAACAGGTCGTGCCGAACAGGAAAAATAGGATGTATGAACGTGCTCATGATTTCGTGACTGAACGTGTGTGCTGGCTTCGGGTCAGTAAGGTCACTGCAAATGAAACTCAGAGGGGCCGGCGGTTAACGTGGCTGGGTCCGTAGGGACACTACGAATTATCACGTCCGCCAGTGGTGACCTTG"));
        // 4. FIB 	Rabbits and Recurrence Relations
        // System.out.println(Rosalind.rabbitRecurrence(5,3));
        // System.out.println(Rosalind.rabbitRecurrence(28,2));
        // System.out.println(Rosalind.rabbitRecurrence(28,3));
        // 5. GC 	Computing GC Content
        // Rosalind.maxGCContent("5_simple.txt");
        // Rosalind.maxGCContent("5.txt");
        // 6. HAMM 	Counting Point Mutations
        //System.out.println(Rosalind.computeHammingDistance("GAGCCTACTAACGGGAT","CATCGTAATGACGGCCT"));
        System.out.println(Rosalind.computeHammingDistance("GTCGCAGGGTATGTAGAGCGTTAAGTGCGGACTATTTCGGCTAGTGTTTTAAGGGCAAGATGCCCTAAAGGTGTAAGTTAGACCCGTGGGAAGAGGTACTCTTACCCGTACCCGCTGGCGCTAGATATTGTGCCGTACACTGGACCTACCACTAGACCCCAATACGGCGAAACGTCTTGTTTACCAGTGTGAAAAGCTGACAATCCACACGATTCCTAAGTCCCGTATGGACGGAACTTGATCTCCCTTGAGGGCATCTGCAACTCATGTTTACTGTTCCGACGGGGATTAAGTGCCGAGAATCTTGACTATGCCCTAGACGAGTCCTCGCGGCTCTTTCGAGCAGTTTTGAGTACGGGCTGCAGACGTTCGAGTGGTTATTCGTGGGCTCGTAAAACTCTACAGACGAACCCCGAGTAAATTTTCTGCCCCGGAGCGACGAGTCTGTGCCGCTACCGGTAAAAAGCAATCTTCACGTAGCGGGATGGCACTAGGTCCGTAGCATTCGTTCGAACTTGACTCGAATGGCAGCGTTGTGCAAACCGAGCTGATGGTAGGAGGTCTCAAGCTTCTTACGCTTATATGACATCCTAACTATGAGACTGTGCTACTTCTTGACATCAGTATTATCATTCTCTGATCAGGTGGATGCCTCGTACTTTAGATGGTCAGCTCATTCAGACTGTTTTTATCTATTTGAGGAGTCCTCGCAGTCGCGTCTCGTTAAACCGCAAGATCGCTGGGACTCTTGCGGCCAATCGTAAGTTCTATCATTAGTAGGCGAATATATCGTCCTCCTAGGCCAGCATCAACTTGTTTCAGTGATGCGGATTCAGCTAAGATCCCACTATGGTAGAGGCTCTATAAAAGATATGGACCATCTTGTCGGCTATGGGCCGCCCCTAGAAAGTGCTATCCTGATGATGTAGCAGTTACACTGGCAAGTTAGGGGTTTTTCCAGTCTTCTCCTGCTCGTTGTGAATGATACTCAA",
                "GTCAAAAAGTCTGTCGAACGTACGAAGCCACCTACTTAGGTGTGTGTTATGACATAAGGGCGTTTACTAAGTGGCAGTTAGGTCGGTTGGGTGGCGTGACCCTTTTCGAAACAGGTCTCGCAACTTTGGGAGACGTCCACTGGTCCTACCAGTTGCAGGTACTACGCCCAAATCTTTCCTTTGCGCCGTCGAGTTGATGACCGACAACGCCCATCTTAAGCTAGGTAGCTAGTAAAAGTGAACCCCGTCGAGTCGAACAGTTATATGTCTATCCTTTTCCGATGGAGACGTAGTCTGAAAGATCTAGGCCATGCCCTGACAGACTCCTCAGGCCGAGTTGAGCCGATTTTGTCTACGCCGTTCCCGTGTTAGGCAGACTTGACATCTTCAAGCAAAATACGAAGTTAAGACTCCGATCGAACAATTGGCATCCGGGTGTAAAGTAAGTGCGGGGCTAAGGCAACTACCATGGTCATGAGGCCATTCCGATGGATGTCTTTCGCATGCGGTTGGAAACGAACCAAATGGGAATGATGTGGTAAGAGTACCGTTGAGAGTCCCTCTCAGACCTATGTCGCTGCTGTTAAATCAATCCACGCGAGATAAGCTGTGTCATGCCTTCGATTATCGCAGTCTCTCAACAACAGGGTATCTCACTTTTGAATTGGTTACATCGTTAAGCGCACTATGCTATGCTTAGGTCGGCCCCGCGGTGGACTCAAAACCACCCAAAGACTCGCCGGGGTTATGCCCGTCATTAGCACTGATGATCATTATTTCTCAAATACGATGCCTTCGTCTGCCGTGATATACAGGTTCAAGGACCGCGGGTGCATGGAAAGTCCCTAACAGCGCAAGCCTCGATAAAAGGTCTCCAACATGTGAACCGCTCTCATCACGTGCAAGAGTTAACTAAGCAGATATGATTCCAAGAACGCTTGTACGCTACCGGTTTCTTTAGTCTTCGATTGCTCAAGAGTAATGGGAAACAG"));
    }
}
