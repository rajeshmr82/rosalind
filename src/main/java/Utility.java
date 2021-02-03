package main.java;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import java.util.regex.Pattern;

public class Utility {
    private static BufferedReader getBufferedReader(String filename) throws FileNotFoundException {
        File inputFile = new File(Objects.requireNonNull(Solution.class.getClassLoader().getResource(filename)).getFile());
        InputStream inputStream = new FileInputStream(inputFile);

        return new BufferedReader(new InputStreamReader(inputStream));
    }

    public static List<String> getInput(String filename){
        BufferedReader reader = null;
        try {
            reader = Utility.getBufferedReader(filename);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        List<String> data = new ArrayList<>();
        String line;
        try {
            while ((line = reader.readLine()) != null){
                data.add(line);
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        return data ;
    }

    public static HashMap<String,String> readFASTAInput(String filename){
        BufferedReader reader = null;
        try {
            reader = Utility.getBufferedReader(filename);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        HashMap<String,String> data = new HashMap<>();
        String line;
        try {
            String id = null;
            while ((line = reader.readLine()) != null){
                var matcher = Pattern.compile(">(Rosalind_\\d+)").matcher(line);
                if (matcher.find()) {
                    id = matcher.group(1);
                    data.put(id, "");
                } else {
                    data.put(id, data.get(id).concat(line));
                }
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        return data ;
    }
}
