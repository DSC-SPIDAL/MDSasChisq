package salsa.mdsaschisq;

import java.io.Serializable;

public class MPIReducePlusIndex implements Serializable {
    public int index;
    public double value;

    public MPIReducePlusIndex(int _index, double _value) {
        index = _index;
        value = _value;
    }
}