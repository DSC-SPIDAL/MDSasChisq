package salsa.common;

public class CommonConstants {
    private static String ERR_WRONG_NUM_OF_BYTES_SKIPPED = "Requested %1$d bytes to skip, but could skip only %2$d bytes";

    public static String errWrongNumOfBytesSkipped(int requestedBytesToSkip, int numSkippedBytes){
        return String.format(ERR_WRONG_NUM_OF_BYTES_SKIPPED, requestedBytesToSkip, numSkippedBytes);
    }


}
