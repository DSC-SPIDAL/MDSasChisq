package salsa.mdsaschisq;

public class FunctionKeyLabelPair implements Comparable<FunctionKeyLabelPair> {
    private double key;
    private int label;


    public FunctionKeyLabelPair() {
    }

    public FunctionKeyLabelPair(double key, int value) {
        this.key = key;
        this.label = value;
    }

    public double getKey() {
        return key;
    }

    public int getLabel() {
        return label;
    }

    public void setKey(double key) {

        this.key = key;
    }

    public void setLabel(int label) {
        this.label = label;
    }

    @Override
    public int compareTo(FunctionKeyLabelPair o) {
        return Double.compare(key, o.key);
    }
}
