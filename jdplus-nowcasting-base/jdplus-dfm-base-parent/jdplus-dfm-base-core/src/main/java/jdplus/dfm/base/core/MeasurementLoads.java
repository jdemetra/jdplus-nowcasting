package jdplus.dfm.base.core;

class MeasurementLoads implements Comparable<MeasurementLoads> {

    public final boolean[] used;

    MeasurementLoads(final boolean[] used) {
        this.used = used;
    }

    @Override
    public int compareTo(MeasurementLoads o) {
        if (used.length < o.used.length) {
            return -1;
        }
        if (used.length > o.used.length) {
            return 1;
        }
        for (int i = 0; i < used.length; ++i) {
            if (!used[i] && o.used[i]) {
                return -1;
            } else if (used[i] && !o.used[i]) {
                return 1;
            }
        }
        return 0;
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append('[');
        if (used.length > 0) {
            builder.append(used[0] ? 1 : 0);
        }
        for (int i = 1; i < used.length; ++i) {
            builder.append(' ').append(used[i] ? 1 : 0);
        }
        builder.append(']');
        return builder.toString();
    }


}
