package salsa.mdsaschisq;

public class ManxcatConstants
{
    static final String PROGRAM_NAME = "MDSasChisq";

    static final char CMD_OPTION_SHORT_C ='c';
    static final String CMD_OPTION_LONG_C = "configFile";
    static final String CMD_OPTION_DESCRIPTION_C = "Configuration file";
    static final char CMD_OPTION_SHORT_N = 'n';
    static final String CMD_OPTION_LONG_N = "nodeCount";
    static final String CMD_OPTION_DESCRIPTION_N = "Node count";
    static final char CMD_OPTION_SHORT_T = 't';
    static final String CMD_OPTION_LONG_T = "threadCount";
    static final String CMD_OPTION_DESCRIPTION_T = "Thread count";

    static final String ERR_PROGRAM_ARGUMENTS_PARSING_FAILED =  "Argument parsing failed!";
    static final String ERR_INVALID_PROGRAM_ARGUMENTS =  "Invalid program arguments!";
    static final String ERR_EMPTY_FILE_NAME = "File name is null or empty!";
    public static final String GRACEFULEND_FILE_NAME = "GracefulEnd.txt";
    public static final String STATUS_FILE_NAME = "Status.txt";

    public static String OutFileNameTag = "manxcat-out";
	public static String ErrFileNameTag = "manxcat-err";
	public static String HtmlFileNameTag = "index";
	public static String StyleFileNameTag = "style";

	public static String Tag = "TAG";
	public static final String TagManxcatConfig = "manxcat_config";
	public static final String TagManxcatName = "manxcat_name";
	public static final String TagManxcatDescription = "manxcat_description";
	public static final String TagManxcatLinks = "manxcat_links";

	public static String SimplePointsPrefix = "SIMPLE";
	public static String ColonPointsSuffix = "Colon";
	public static String GroupPointsSuffix = "Group";
}