package salsa.configuration.sections;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

public class MDSasChisqSection {
    public MDSasChisqSection(String configurationFilePath) {
        Properties p = new Properties();
        try {
            p.load(new FileInputStream(configurationFilePath));
        }catch (IOException e){
            e.printStackTrace();
        }
    }


}
