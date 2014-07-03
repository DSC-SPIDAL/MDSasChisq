package salsa.mdsaschisq;

import com.google.common.base.Strings;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

public class MetaDataIO
{
	public static void WriteResults_Cluster(String ControlFileName, ArrayList<String> linesToOutput)
	{
        if (!Strings.isNullOrEmpty(ControlFileName)) {
            try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(Paths.get(ControlFileName),
                    Charset.defaultCharset()),true)){
                writer.println();
                for (String line : linesToOutput) {
                    writer.println(line);
                }
                writer.close();
            } catch (IOException e) {
                SALSAUtility.printAndThrowRuntimeException("Failed writing data " + ControlFileName + " "  + e);
            }
        }
	}
}