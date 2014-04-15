package salsa.configuration;

import salsa.configuration.sections.MDSasChisqSection;

public class ConfigurationMgr {
    private String configurationFilePath;
    public MDSasChisqSection mdSasChisqSection;

    public ConfigurationMgr(String configurationFilePath) {
        this.configurationFilePath = configurationFilePath;
        mdSasChisqSection = new MDSasChisqSection(configurationFilePath);


    }

    public static ConfigurationMgr LoadConfiguration(String configurationFilePath){
        // TODO - Fix configuration management
        return new ConfigurationMgr(configurationFilePath);
    }
}
}
