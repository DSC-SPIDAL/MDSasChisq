package salsa.configuration.sections;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.regex.Pattern;

public class ManxcatSection {

    public static void main(String[] args) {
//        generatePropertyRead();
        generatePropertyWrite();
    }

    private static void generatePropertyWrite() {
        String file = "/home/saliya/Desktop/manxcatsection5";
        try (BufferedReader reader = Files.newBufferedReader(Paths.get(file))) {
            Pattern pat = Pattern.compile("[ ]");
            String line = null;
            while ((line = reader.readLine()) != null){
                String [] splits = pat.split(line.trim());
                String var = splits[0].trim();
                String varName = splits[2].trim();
                boolean isOtherThanString = !varName.startsWith("p.");
                int start = varName.indexOf('"');
                int end = varName.indexOf('"',start+1);
                varName = varName.substring(start,end);
                String prefix = isOtherThanString ? "String.valueOf(" : "";
                String suffix = isOtherThanString ? ")" : "";
                System.out.println("p.setProperty(" + varName + "\"," + prefix + var + suffix + ");");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void generatePropertyRead() {
        String file = "/home/saliya/Desktop/manxcatsection4";
        try (BufferedReader reader = Files.newBufferedReader(Paths.get(file))) {
            Pattern pat = Pattern.compile("[ ]");
            String line = null;
            while ((line = reader.readLine()) != null){
                String [] splits = pat.split(line.trim());
                String type = splits[1];
                String var = splits[2];
                String prefix = "int".equals(type) ? "Integer.parseInt(" : ("double".equals(type) ? "Double.parseDouble(" : ("float".equals(type) ? "Float.parseFloat(" : ("boolean".equals(type) ? "Boolean.parseBoolean(" : "")));
                String suffix = !"String".equals(type) ? ")" : "";
                String def = "int".equals(type) ? "0" : ("double".equals(type) ? "0.0" : ("float".equals(type) ? "0.0f" : ("boolean".equals(type) ? "false" : "")));
                System.out.println(var + " = " + prefix + "p.getProperty(\"" + var + "\",\"" + def + "\")" + suffix + ";");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
