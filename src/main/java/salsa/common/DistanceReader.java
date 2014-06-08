package salsa.common;

import com.google.common.io.LittleEndianDataInputStream;

import java.io.*;
import java.nio.ByteOrder;

public class DistanceReader {
    public static short [][] readRowRangeAsShort(String fileName, int startRow, int endRowInclusive, int globalColCount, ByteOrder endianness) throws IOException {
        int numRows = (endRowInclusive - startRow) + 1;
        FileInputStream fis = new FileInputStream(fileName);
        DataInput di = endianness == ByteOrder.BIG_ENDIAN ? new DataInputStream(fis) : new LittleEndianDataInputStream(fis);

        int numBytesToSkip = startRow * globalColCount * Short.BYTES;
        int skippedBytes = di.skipBytes(numBytesToSkip);
        if (skippedBytes != numBytesToSkip) throw new IOException(CommonConstants.errWrongNumOfBytesSkipped(numBytesToSkip, skippedBytes));

        short[][] buffer = new short[numRows][];
        for (int i = 0; i < numRows; ++i){
            buffer[i] = new short[globalColCount];
            for (int j = 0; j < globalColCount; ++j){
                buffer[i][j] = di.readShort();
            }
        }
        return  buffer;
    }
}
