package main.java;

import main.java.biostronghold.Rosalind;

public class Solution {
    public static void main(String[] args) {
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
        //System.out.println(Rosalind.computeHammingDistance("GTCGCAGGGTATGTAGAGCGTTAAGTGCGGACTATTTCGGCTAGTGTTTTAAGGGCAAGATGCCCTAAAGGTGTAAGTTAGACCCGTGGGAAGAGGTACTCTTACCCGTACCCGCTGGCGCTAGATATTGTGCCGTACACTGGACCTACCACTAGACCCCAATACGGCGAAACGTCTTGTTTACCAGTGTGAAAAGCTGACAATCCACACGATTCCTAAGTCCCGTATGGACGGAACTTGATCTCCCTTGAGGGCATCTGCAACTCATGTTTACTGTTCCGACGGGGATTAAGTGCCGAGAATCTTGACTATGCCCTAGACGAGTCCTCGCGGCTCTTTCGAGCAGTTTTGAGTACGGGCTGCAGACGTTCGAGTGGTTATTCGTGGGCTCGTAAAACTCTACAGACGAACCCCGAGTAAATTTTCTGCCCCGGAGCGACGAGTCTGTGCCGCTACCGGTAAAAAGCAATCTTCACGTAGCGGGATGGCACTAGGTCCGTAGCATTCGTTCGAACTTGACTCGAATGGCAGCGTTGTGCAAACCGAGCTGATGGTAGGAGGTCTCAAGCTTCTTACGCTTATATGACATCCTAACTATGAGACTGTGCTACTTCTTGACATCAGTATTATCATTCTCTGATCAGGTGGATGCCTCGTACTTTAGATGGTCAGCTCATTCAGACTGTTTTTATCTATTTGAGGAGTCCTCGCAGTCGCGTCTCGTTAAACCGCAAGATCGCTGGGACTCTTGCGGCCAATCGTAAGTTCTATCATTAGTAGGCGAATATATCGTCCTCCTAGGCCAGCATCAACTTGTTTCAGTGATGCGGATTCAGCTAAGATCCCACTATGGTAGAGGCTCTATAAAAGATATGGACCATCTTGTCGGCTATGGGCCGCCCCTAGAAAGTGCTATCCTGATGATGTAGCAGTTACACTGGCAAGTTAGGGGTTTTTCCAGTCTTCTCCTGCTCGTTGTGAATGATACTCAA",
        //    "GTCAAAAAGTCTGTCGAACGTACGAAGCCACCTACTTAGGTGTGTGTTATGACATAAGGGCGTTTACTAAGTGGCAGTTAGGTCGGTTGGGTGGCGTGACCCTTTTCGAAACAGGTCTCGCAACTTTGGGAGACGTCCACTGGTCCTACCAGTTGCAGGTACTACGCCCAAATCTTTCCTTTGCGCCGTCGAGTTGATGACCGACAACGCCCATCTTAAGCTAGGTAGCTAGTAAAAGTGAACCCCGTCGAGTCGAACAGTTATATGTCTATCCTTTTCCGATGGAGACGTAGTCTGAAAGATCTAGGCCATGCCCTGACAGACTCCTCAGGCCGAGTTGAGCCGATTTTGTCTACGCCGTTCCCGTGTTAGGCAGACTTGACATCTTCAAGCAAAATACGAAGTTAAGACTCCGATCGAACAATTGGCATCCGGGTGTAAAGTAAGTGCGGGGCTAAGGCAACTACCATGGTCATGAGGCCATTCCGATGGATGTCTTTCGCATGCGGTTGGAAACGAACCAAATGGGAATGATGTGGTAAGAGTACCGTTGAGAGTCCCTCTCAGACCTATGTCGCTGCTGTTAAATCAATCCACGCGAGATAAGCTGTGTCATGCCTTCGATTATCGCAGTCTCTCAACAACAGGGTATCTCACTTTTGAATTGGTTACATCGTTAAGCGCACTATGCTATGCTTAGGTCGGCCCCGCGGTGGACTCAAAACCACCCAAAGACTCGCCGGGGTTATGCCCGTCATTAGCACTGATGATCATTATTTCTCAAATACGATGCCTTCGTCTGCCGTGATATACAGGTTCAAGGACCGCGGGTGCATGGAAAGTCCCTAACAGCGCAAGCCTCGATAAAAGGTCTCCAACATGTGAACCGCTCTCATCACGTGCAAGAGTTAACTAAGCAGATATGATTCCAAGAACGCTTGTACGCTACCGGTTTCTTTAGTCTTCGATTGCTCAAGAGTAATGGGAAACAG"));
        // 7. IPRB 	Mendel's First Law
        // System.out.println(Rosalind.probabilityOfDominant(2,2,2));
        // System.out.println(Rosalind.probabilityOfDominant(20,17,17));
        // 8. PROT 	Translating RNA into Protein
        // System.out.println(Rosalind.translateRNA("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"));
        // System.out.println(Rosalind.translateRNA("AUGAUGUCGCAUAUUCUAAGGUGGCAUGCUCAUAAAUUGCCCAUUGUAUGGAAAGUUAGAAUUGCGAUUGACGCGCUUCCCAAAACUCGAGCACGUCAGUCGGAAUCUGAUCAGGCGGUCCAUCUUGUGGGAGGGAACUGCGAAUUCCGAAAUUCUUUAAUGGACGAUGUCACUAUCCGAGGGCUUCAAUCCCUGGGGCAAAAGCGGGGGACACCGAGUCUGCUUAAAUUAGAGCGAGAUGAACCCGACCCCAUGCGUCGACUCGUGUCAACUGAGUCUAAAAGGCCAUCUAGCUGCUAUCAAUAUCUUCUCCGCGUCCAGUGGAUGACCACCGGCUUUGUUACUUGCUGCUUCGGUGGCGCCCACAAAUGCGUUCGCUAUUGUUCAUUGUGCGGGGUUUUGGGGACGUCUAAAUCCACAACCGAAUUCAAGUCUACGAGAACACCACAGCAGUUGGCAGCGGCGUUCUCUUCUGGGCUUCCCUUUGGGAGGCAGCGACCACCCCCGUGCCGGUGUCAAGUUCGUGCGUGUUUUGGCGGCAGUCUUUCCAGCCGCCCCGUUGGGUUUCAAAAAUUGGGGGCUAAUUACUCUCUCGUGUGUUUCUCUACUAUGUCCUCCGCAGGGACAGUCGCCCUGGUCUCAUCGCUCACAAUCUUUGCGCACUCAGAGGUUGGCUACAUGGCUGCCUAUUCAUCAGCACACCGCAAAAACUUUGUAGGGGUGUCCCGCUUCCAGUCGGAGGGUAUAGGUCGGGGUUGUCGUCAUUACCUUUAUUUCCACCGUCAUUAUCAAAUGGCAUCUUACGUGGAGAUUCCAACCACCGAGUUUCGACCAACCAAAUUACCCGGUUGUAUCCUGACCGACGCUGAUUUGGCCAUCUGGGCGUGCGGUGCUUUCCCCCGUAACUAUCUGAUUAUCAUCGACUGUGUGUCUGCUAAUCAUUGUGCUAAGAAAAGUCCAUCAUGCCCAUUAAAGGUACAUCCUGGGCGGUACACAGGCAGGCCUUCCGCUACGGCCGAGCUUCAUGCCUGGUCACAGAUGAAGCGGCUCGCAACGCCCUUGGGACUAGCGGACGGUACAAAGUCUCAGACGCCGCCGGGGCAAUAUUACUUGCGAUUACAGUUUGAGAUCGGUGGCGAUCCUCCGAGUUCGAUAGCCCCUGUAACACCUUGCAAACACCGGCGACCACACGGUACGUUUUCGGAAGAGCCAUUCGCGUAUAUGGCGAUCAAUAAUGCCGCCAAAGCUCUAAGUGGCCUUCUUACUCAGCGCAAACCAAGGGUUCCUACGCCAUCUGUGUGGAUGCAGAGUGAGUUGCAUGUUCUAAUCAAGACACUUCGGACUGUGAUAGUGCUGACCAUCUCCUUUUGCUACCUUUAUAUAAGAAUAAUGCCAGUAGCGGGGCAAGUCCAGCCCAUGUAUCUGGCAUCUGAAAAGACUGCGAUCUUAUGCCCAAGAUCAAAGAGCAGAGAUGAGUUUUUAACGAACGCGACGUUGAGGGUAUUAUAUAAGACGAGCUCGAUUCGCUUACGCAAGAUUCUUCACGGUUUAAGUCUGGGCUUCCUUUCCGUGUCCUGUGCACGUAGACACGUGUUCGUCCCGCACUUGGUGGCGAUAUCCAGUCCUUCGCCUUCGAUACGGUGGACCCGGCAAUACGUACAGCGCUACCUCUUUCGACAGUACAGGGUGGGAUCCCAUUGCCAUGAUCUAUAUUGCCGGGUAUUAAUCCGCUCCUUCGGAGACACGGGUAUAAGUGCGAUUCAAAACGCCGUGGUGCUGAUACGGGCGAAGGGUGACUCGAUGAUAAAUACGGCAAUCCAACUCCGCGUAAGUGGGGGCCCAUCCACGCUUAAACGAUCUCUCAGGGACGCAGAGGCUCCAAAGUCUUCUUAUCAUGUAUCUAUAAUUAAUGCCCGAGCGGUAGAUAUAGGUCUUAUGCAAUAUACAGGCAAGCUAUGCAAGCAGCAGUCAAUCGCCCGCAGAAAUACCAACCCCCCUUCCCGUAGGAAGAUUCAUGCACCAUGCGGACUCCGUCUCGCACCCCUCCCAUCCAAUUCGAGCGAGAAGGUCGUCGUGAACGUUUCCUUUCUUUCCCUGCAUGGGAGACCAAGACGUCAGGUGCCUACCAUUAACACUAGGCCUUGCGCUCGCAACUUUCGUGUGGCUGUAUGGGUACCUUCUGGCUCACUACCUCUGUCGGGGACAUUCCAUCUUCAAAGAUGUUAUAUACCGGAAGUAGUAACCUGCAUUGUUCGGAAGCAAUCGAACAGUUCCGAUGGUACAUGGUACCCGCAAAUCACGCGUAUGAGCGUAAUGUACAACGACAGUUGUUUGGUGACCAGGCAUACAGGAUCAAAUAAUCAAAUUGUUGAUGGGAAGAGUGGGCGAACUGCGGAGGCCCCAUUCGCUACAAGGGUUCUUUCGAGCUGGCCUGAAAAGGCUGUUCUCACCACGUUUUGGUGUCAUAAGUCCUCAGCAGAUAAGCUAGGCUCAGCGCUCGUCGGAAAUGCGAAUUCUGAUUGCCAUUCUCCGUACACCGCUCCUAAGUCUUCCUUUGUGCCCGAGAACCAAUACGCAGCCAAGACGACGCGAUCACAGCAUCUACUCCUACCUACACCAGAUUGUAUAUCAUACAGAUACGGAGCUAAUGCAACGGAGGCCAAGGGUCACAAUAUUCGGUCAUCGCGGAAAACACCAAGGACCCGGGAUGCCACAAUACAGAGACAUGUAACCGGUCGCGGACUCCGCGAGAUGCGUGGCGAAUCGGCUACGCAUCGGCGUAUAAAUAUACAAGCCCCGAGGCACUCUCGAGCAAUAUCAGGGAGCGGUCCGCCCUCCAUUAGUUCCAAAGGCUGCUACACACUUAGAGUUCUGAUAUUUGAUAUCGCCUGGAUAUCCGGCAUCCGCAGCCAAGUUAAACGCGGUAACGGUUGUCACGAAGCAGCCGCAAUACACCCACGCAUCUUCACUCGCGGUCUGCAACUUUGUGGUUUCUAUUCGUAUAUUCAAGGACUAGGAGAAGUAGACUGUACAUCUACCAACAUAGCAUCCCUCCACCCUGGGAACGAGCGUUCUAGACGAAGCAGACCUGCCCGCGUAGAUGUUUACGAGUCUCCCUUCCCCUUUCACAGUGUCGGGGAUUCGAGACCACCUCGGUUGGCAAUGGCAUGCCCGAUUGGAGAGCACUGUAUGAUGCUGCGUACUCUUUGCCGGUCACUCGAUCCGUAUCUUGAUUUGCAACCUGCGCUAUGCAUACAAACGAAAUGUCGACCGAGGAACCGGAAAACACGUCUGAAACCUACCUUCGGUCCGAAUGUACCUGAGAGCGCAGUAGUUGGGAACGAACCCCGGGGCCUCGCUUCCCUACGCAUUAAUAAACGAUUGUCUCAUCUGAGGAGUAUCUCGUAUCCUACGCCUGAACGAUUAUCUCAACCGUAUCGACCCGUCAGGGGAAGGAAGAAUGGUUCAACUCUAUGGAUUGAUUAUAGAGAUUUGUCAGUGGUAGUGGAUGAUACGGGACGAGUGUUUGGGCAUAUCGUACGAAUUCAAUUGUCUGAUAAUAUGGUGCUUCAGCGUGGCAAGACCCAUCGAGUCGCACAGGGAUUCCAUUGGCAGAAGCGAUACAGGCGCGCAUUCACCUAUAACGUUCAUAUUAUCAUAGGCGCGAGGGCCUUUGUAGUGUUCGCCCCGGACUGCGGGUUCGUUCACAUUUGCAUCAAGGCUUUCACGACUUUAACGUACCUAUGCCGGGUUGUGAUAGUUUGCAGCACUUUUACGCACUGUUAUGGGCUCUACAUGACGGUGGCAGGGUAUGUAGUUAUGAUGACAGCCUUCACUCUCGCCGGCAUGCGGGACGUCUUUGAAAGGCUCAUUGCAUUAAAUCCUUCCCCCGUAUCGCGCAACUAUGUGCUUCAAGGUAUAAUAAUUCGCGUGGAGCCGGGAAGGUCGGACGUAUGGUGGCGAUAUCCAUCAGACUUCAUGCGGUCCAGUCGCCCGAAGCAGUUCGAGCAUGAGCUCUAUUUGGCUACUGCGACGCCUCCAACGCCUUUAAAACGGUCCGACCCGUGGACUGAUUGCGGUCAUGGGCGCUUGCCGGGGAAAGUGGAGUACCCAACCGCCCCGAGCCCUCACGCCUUCGGGUUGCAACGCACGGAUUUCUCUGCUGUGAAAAGCACUGCCCCGGCUAUCCACUUGCGAUCAUUUAUUCGGAGCGAACAGUUAUCCAUCGCGCUACGUCGAGGACUUGGGAUUCGUUCAUGUCAGGGACGGACCCGACAUAAGAUGUUCUACGAAGGCAGUGGAUCAUCUAUCACGCAUAGAUUGUUGCCGUUUAACGAUACGCCUCGAGGCAACUUGAGUCAUGGCGUACCACCCAGCCGCCUUGCAAGGAACUGGGCUAUGCUUAGGUGGAGCUCGUGGGUCAUAGCUCGCUCGAUAACGCCCAGAUUGUACCACAUCCAAGGAACUACGUUUUCUUCGGCUCGUUUGGAAUCCGCACUCUUCCAUGUAGGUGAAGUGCGCGGGGCCAGCUCAGCGCGCCCGAUUAGAAGGAGACCUCCGAGGCUCUAUGCGCUAAGUAAAAUUCUCCUUUUACCCUGCAAUCUAACAGCUUAUUAUUCCUCAGAGGGAUUUUCACUCGUGCCGGGCUUCUUCGUCAGCCGUGCAGUCCCUAACCAAUCGCAACGCUUAAAACAGCUAUCUCGCCUUGCUAAGUCCCGGAGUUUGUGUAGUCAAAAUCAGUUCCUGGCGUGUCCUGACGUAAUUUUAAUGGGUCAAGUUCCCGUAGAUGCCAUGCACCAACCGCCAUGCUUUGUCACCUACAUAGGCAGUGCCGAAAGCAAGCGAUUUAGAAUGCAAUAUCAAUCCGGCGGCCGAGAAAUCCUUUAUCAAGGAGCACAAGGUGCGCGUAGACCGCGUAUAUGUCAUAAGCUAACGAUAGCCAGCGACCAGACCAAGGCAGCGGUAGAGCCUGCCAGUAAGGCAUGCCGUAUGAUGCUAUCCUCAGGCCCGGCGUGGUUACACCUACUGCUGGGGCCCAAGCGACCACUUGUGUUAGUGGGACUGCUCCGGUUAUCUUCGUGUGUUCCAUGUUGCCGUUCGUGGCGGCUUGGGCCACGUGAGAAGCCUAUCAGCUGUGUUGGAAUUUUUAAACACCGCUCCGGCGACCCACUACACUCCGGCGCGGUCGGCUUGGAGCACACGCUCACAACACUAUCUCCUUUCGGGCAAAAUAAGCCACGGACCCAGGCAUACGUCAGAAAUCGUUCCCUGUGUGCGAAACGAUCCUAUCCACCCCUGUCUCCAGGAUUGUGUGUCGAAAAUAUGGUAACACGCUUUCUUACCCGCCUAAACAGGACUAUCUACUUCCCGCGCAGUGUUAGACUCUUUCAGUCGACAGGGGGGGCGUCCAUACGCUUGGUUUUUAGCCAGCUCCUGGGCGUGAAAAUACACAUGAUAUUAGAACCGUGUGUUUUGCUGAAGCAAGAUCCGGGCUGUGGAAGAAUACUCCCAAAUCUUUGGGACGCUUCGAGAUGGUUUCAUCUCCAGAUUCGAUGGGUCGGUAAACGUCGGAAGGAAUCUCGAGAUGUUCCCGAGACCGAGCGCUCUCUGCACGUUACUGUCCGCAGUGUAGACAUGGGUGACACCGGAUCUAGCUAUUGCCGUGAACUUGGACACAAGUUAUUCCAUAUCGCGCAAUCGCUAGUCUCGGGGGGUCUACACGAGUUAAUCUGCAACGAUAAUGAAAAUAGUGUAAGAAUCUCACAGGAUCUGCUUUGUGGCUCAUUUGAUACCCCCCGUUCAUACAAUGUUAUACUAGUAGGUUUCUGUACUGUCCGGUUUUCUGUACGAAUGGGCUUACACCCGAAAGAAGCACGAUCGUGGUAUCUCCCCCUAAAGCCACACGAUAUACCCUCAUACUUCACUGAGGGUGCGAAAGGACAGACUGGCCCGAAUACAAGUGCCGCGUAUUUGUUCCGGGUCAAUUCCGUUCGUGCUCGCUACCAAAUUCUGCAUGGCACGGGCACCUGGGGCAGGCGUGGCCCGUGUAUUCUGAUGCACCAACUAAGCAAAAGCAGGGCCGGUCGAGACGGUAUCACAACGAGGAAUCGUGUCUCCUCGUCCACCGAUUGUGGCAAGCAGGAAGGAUGCCCUGCUAAAAUACAUGAUGCCCACUUGGCUUGCCACAAACCCAAACAACACAUUUCUUGUCCUUCACACAGCCACUCGGAUCACAGGCUGAUUUUUCAAACCACAUAUGCUUCCACUGCAUUACCACGCCGUUCGUAUCCCUAUAACCACUUGUCUGAGCUAAUUAUCACCCUUGGCGGUUAUUCAAUAGCUCUUAAUUGUGUAUCAGGAAUCCCCACGCCGCACAACGUUACUGUACUUAUUAUGAGUUCGAGCUUACAAAAUGGAGGGUGCGGAGUUAAAGCCGAUCAUCCGUUUAACGCACACCGGAAGGUUAUUGGUCAGAUUUCUAACGUGAUUAACGAAUGUAUAGCUGCGGAUUUACGCUCAGCGGGGGAUGGCCCGAGGUUCAGAGAGAAUGAGGUACGUGGCAGAGACAGCGAACUCGGUGUACGCACACAGCAUCGGUUCAGGCUUCGAACACCCUUUGAGCUCGCGCUCACACUCAGCCCCAAUCAAAUGCGCGUCAAGAGGGACCUACUCACACGUCAUUACGACCAGUGCUACUUUAGUUCCUGCAUCGAUGGAGUUAGCCCGCUUAAGACGGAGCGACCUCCUAAUCACCUAGAUGGGAAAUUUCCCUUCUCGGAUUCGCUACCACUCGUAUCGGUGUAUCUUGAUAGCAGCGCACCGCUCUACGGGUUCGUUUACGCCAGUCGUAAGGACUCCAGGGAGCCGAGCUGCCCGUGGCAAGAGAAAUCCAUAAAGAGAUAUACUGGUAUAGACAUGGACGCUGCCAGUAUAACACGCUGCGUAAUAAACCGAAGCACCAGAUGUGCUAUACUCGGUUCUCCUGUCACCGUAACUACUGCUCGGGAUACCGCGCCCCGUAACGCUGAUUACGCGGCAAUUAUCGCAGUGCUCACCCGUCCACAUACAAUGGACCUACGACUACGGAUUUUAGCGGGAGGCGUAUAUACGGGCAGAGCCCCAAACGAGGUCAUGGUACGUCGGGGAAAUAACUCUAUAUCUUACCACCUGUUUGGGCCUCACUUACUUCUCUGUUUUAGUAUCUUCGUUGCCUGGGUCAGAAUUUCGGGGCCUCGCUGGCCGGCCUGUCCAUCGCUGGUACGCCUGGGACACAAAUUACACACGAUUCGCCCGUAUAAGCGUGACCUGUUCGCGGCGCCAGCGAUUUUGAUUCCCCAUCCGAUAAGAGGUGUGUGCGCGAGUAGCUUAUCCUGUAAGGCUAAAGCUGAAGAUCGCGUGUUGAAUGUUCACACCGAUAAACAAGACAGUUAUUCGGAAAUUUCUGUCAAUCCGAGAACUCAGUUUCGCUACAGAGACUCCGUUGGAUCACCCCUCGAAUACAUAGUCCGGAUUUCUCCCCGCGCUUAUACUUCAUUUGUCAAUCUGUGGUGUCGUCUCCAUGGACUUAGAGCCCCACCUAGCCUCAUAGUACCCGUCUUUCCCUGUGUAUUCCUCAACACCUUGAUCUACUCAGUUAUCCAAAGGUUCCAGAAUAGUAGAUGCGUUUUCCGCCCAGCCAUAGAAGUACGUUUGCGACCUUUCGCCCGCACGGCCUGGCGCAACCAAUGCCCCGUGCCUCUGUUACCGAUGUUAACGCUUAAUUAUUGGCAAAUUAUCCUUCAUACGAGGCCUGAUACAAGCUCAGUUACUCAUGACAAACCAAAUCUAAACUGUUGCAUAGGGUAUCAAGUAAUAUGGGACCGCCCAAUAGGAGUUGGGCAACAGCUGCAAUUCCAGACCUUGCUAGUAGAAUCAAGACCUUGCUGCAUAUUGUCUGUCACCCCAUUCUGUCGUUUCACCAGCGAUUCGAGCGCAAUGUCCAUUUGUGCCCCACGAAGUAGGUGUCCCCAAGCGCCCAUUGAGAACAUGGUGUUGCUGCUAAGAAUCACAACAGACCUUUACCGCCGGAAUAGGAUAAAGUCUGCCGGGUACGCUAACGCCUUUAGUGAUGUGCAAGAGGUUGUCCAACCCGACACUGGCAGGUGUUUGCGAGAUUUGAGCGACCAGCCUGAGGGCUGUCUAGUGUCGGGCAGUUUCGUGUACCCGUCAGAGCUUCAUGCGUCACAAAUGUAUCAUCAGAUUUCAGAACUCGGCCUCCCGCCCUUGUAUGACCCCCUGGGGAGCCCAACAUAUGUCAUCCCCAAACAGGCCGUUACAACAUAUUCAUUAGGUAAAGGAUUUUCGUCGCCAAUACGCGUGACGGCCUGGCGAGUUCAUAGGCCCGUUACAUGUCGUGGCGACCCAAGUCAUCCCCUAGGUGUUGGCUUCCCGUAUGACGAACCGGGAGCCUAUAAGAACGGGUUAACGCUCAUCGAAAGGUGGCAUUUCGUAUACCCAUUACCACCUGAUUUAAAUGCAGUCAGAACGUCAGAUCGAAAAGCACAGUCGACACCUAUACUCCCCCCUAAGUCUGCUAAUAGAGAGCACAUCAGUGGAGGAAUCAAGACGGGAUGCUCGCCGACAACUAGUCUACAACGGGGCUCUAUAGAGCUACCCUGCUUGAAGCGGUGCUCGCGCACUAUAUCUGGCCUACAGUUCAAAAGAUUUAGGUCGGGAGACAUUUGUAGCUCAAGAGUAAACGUAAGGAAAGCGGGGUCCGAUUUUCAUGCAUGGGUGCGCUUUACGGUGGCGCGGGAACAUGCUGGGGAGGUUCUCAUACAACAGACGGUUUUGUGGAGAUUCACGCAAACUGCAUCGGGGGAAGAUUACUGCGAAUUUCUACCCGAGACAUAUCCGGAGGACCAAAGUAAUAGAAGAAGCAAUCUGCGACGAUGUAGCUGUACCGCUCAUAGGUUACAGUCUCUUAGCUUUGUAGGACCGCAAGAUCACGUUGGCCGGAGCCCGCGUGCUGGUUACGUCGUAAUACCCCAUGGAGUCGGGGUUAGAUAUGGGAGACGGGCGACGUUACCAUCAUUGUGCGCCUGUGAGGUCCAUCGAUCUAUAUCAGACACUCGACCGGCCCAUAUAGAAGACGUGUAUAUAAAAAACAAUUCAUGCUAUUGUUGGCCCUCGGCUUUAACCCCUCACAAGGCAUCCAGAUCAAGUUCGCUGUUACAAACUCCGACCGGUUUACAUCCAAAUUACUCAUGGCAGGGUCACCUCGGUAUGCAUCUGCAUUUCGGUUCCGCGGACCGCGAUUUCACUAUGUACGUCCUACGCGGUGUUCUUAGCCAUCACGUUCAGGGCCUAGCAACAUUAGUCGUUAUAAUGGUCCAAUAUGUAUAUUUUCAAUUGGCCCCGCACAACGGACUCAAUGAGCUGAAAACACUUAGUGAGUGA"));
        // 9. SUBS 	Finding a Motif in DNA
        // System.out.println(Rosalind.findAMotif("GATATATGCATATACTT","ATAT"));
        // System.out.println(Rosalind.findAMotif("GCCCCCAAGCCCCCAACTCCCCCAAATTTTACCCCCAACAGACCATAAGCCCCCAACTCGCCCCCAATAGTGAGCCCCCAACAGATCCCCCAAACCCCCAAGCACCTCCCCCCAAGTCCCCCAAACCCCCAAACCCCCAACCCCCCAACCCCCAACCCCCAAAACCCCCAAGCCCCCAACCCCCAATCCCCCAAGGGCTCCCCCAAATCCCCCAAGACCCCCCAAGCCCCCAATCCCCCAAAGTACCCCCAAGACCCCCAACCCCCAAGTGCGGATCTGGCTACCCCCAAGCTACCCCCAAGTGGCCCCCAACCCCCAACCCCCCAACCCCCCCAAACCCCCAACCCCCAACCCCCAACCCCCAATTCTCTTCCCCCAATTCCAAATCCCCCAATCCCCCAAGTCCCCCCAACAGCCCCCAAGCCCCCAACTGCCCCCAAGCCTCCCCCAACCCCCAAATCCCCCCAATCCCCCAAGGTCCCCCCAACCCCCAAGACCCCCAAGATAACACCGTAGCCCCCAACTTGCCCCCAATTAGGACCCCCAACCCCCCCAACGCCCCCAACCCCCAACCGCCCCCAAAGTTGGTGGTCCCCCAATTGAATGAAGCTACCCCCAAAGTAACTAGACCCCCAACACTACCCCCAAGCCCCCAACCCCCAACCCCCAACCCCCAACACCCCCAAAAGGGCCCCCAAAGTATCCCCCAACTCCCCCCAACCCCCAAGCCCCCAATATCCCCCAATCCCCCAAGGGCCCCCAAGAGATTCCCCCAACCCCCAACCCCCAAATCCCCCAACCTTTCCCCCAATCCCCCCAA","CCCCCAACC"));
        // 10. CONS 	Consensus and Profile
        // Rosalind.getConsensusString("10_simple.txt");
        // Rosalind.getConsensusString("10.txt");
        // 11. FIBD 	Mortal Fibonacci Rabbits
        // System.out.println(Rosalind.mortalRabbits(6,3));
        // System.out.println(Rosalind.mortalRabbits(83,16));
        // 12. GRPH 	Overlap Graphs
        // System.out.println(Rosalind.overlapGraph("12_simple.txt"));
        // System.out.println(Rosalind.overlapGraph("12.txt"));
        // 13. IEV 	Calculating Expected Offspring
        // Rosalind.calculateExpectedOffspring("1 0 0 1 0 1");
        // Rosalind.calculateExpectedOffspring("18960 17972 16529 17477 19039 19698");
        // 14. LCSM 	Finding a Shared Motif
        Rosalind.findLargestCommonSubString("14_simple.txt");
        Rosalind.findLargestCommonSubString("14.txt");
    }
}
