/*
 * Copyright 2023 National Bank of Belgium
 * 
 * Licensed under the EUPL, Version 1.2 or – as soon they will be approved 
 * by the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 * 
 * https://joinup.ec.europa.eu/software/page/eupl
 * 
 * Unless required by applicable law or agreed to in writing, software 
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and 
 * limitations under the Licence.
 */
package demetra.dfm.timeseries;

import demetra.timeseries.TsDomain;
import demetra.timeseries.TsPeriod;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 *
 * @author Jean Palate
 */
public class TsInformationUpdates {

    /**
     *
     */
    public static class Update {

//        Update(final TsPeriod period, final int series) {
//            this.period = period;
//            this.series = series;
//        }
        Update(final TsPeriod period, final int series) {
            this.period = period;
            this.series = series;
        }

        /**
         *
         * @return
         */
        public double getObservation() {
            return y;
        }

        /**
         *
         * @return
         */
        public double getForecast() {
            return fy;
        }

        /**
         *
         * @return
         */
        public double getNews() {
            return y - fy;
        }

        /**
         *
         */
        public final TsPeriod period;
        /**
         *
         */
        public final int series;

        public double y, fy;

        @Override
        public String toString() {
            StringBuilder builder = new StringBuilder();
            builder.append("var:").append(series).append('\t').append(period)
                    .append('\t').append(y).append('\t').append(fy);
            return builder.toString();
        }
    }

    private final List<Update> news = new ArrayList<>();
    private final List<Update> revisions = new ArrayList<>();

    TsInformationUpdates() {
    }

    /**
     *
     * @param p
     * @param series
     */
    public void add(TsPeriod p, int series) {
        news.add(new Update(p, series));
    }

    public void addRevision(TsPeriod p, int series) {
        revisions.add(new Update(p, series));
    }

    /**
     *
     * @param freq
     * @return
     */
    public TsPeriod firstUpdate(int freq) {
        TsPeriod first = null;
        for (Update update : news) {
            if (first == null) {
                first = TsUtility.lastPeriod(update.period, freq);
            } else {
                TsPeriod cur = TsUtility.lastPeriod(update.period, freq);
                if (cur.isBefore(first)) {
                    first = cur;
                }
            }
        }
        return first;
    }

    /**
     *
     * @param freq
     * @return
     */
    public TsPeriod lastUpdate(int freq) {
        TsPeriod last = null;
        for (Update update : news) {
            if (last == null) {
                last = TsUtility.lastPeriod(update.period,freq);
            } else {
                TsPeriod cur = TsUtility.lastPeriod(update.period,freq);
                if (cur.isAfter(last)) {
                    last = cur;
                }
            }
        }
        return last;
    }

    /**
     *
     * @param freq
     * @param updates
     * @return
     */
    public TsDomain updatesDomain(int freq, List<Update> updates) {
        TsPeriod first = null;
        TsPeriod last = null;
        for (Update update : updates) {
            if (first == null) {
                first = TsUtility.lastPeriod(update.period,freq);
                last = first;
            } else {
                TsPeriod cur = TsUtility.lastPeriod(update.period,freq);
                if (cur.isBefore(first)) {
                    first = cur;
                }
                if (cur.isAfter(last)) {
                    last = cur;
                }
            }
        }
        if (first == null || last == null) {
            return null;
        }
        return TsDomain.of(first, first.until(last) + 1);
    }

    /**
     *
     * @return
     */
    public List<Update> news() {
        return Collections.unmodifiableList(news);
    }

    public List<Update> revisions() {
        return Collections.unmodifiableList(revisions);
    }
}
