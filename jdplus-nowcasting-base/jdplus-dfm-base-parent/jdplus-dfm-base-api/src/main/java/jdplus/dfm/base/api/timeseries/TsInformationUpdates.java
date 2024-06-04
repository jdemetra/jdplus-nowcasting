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
package jdplus.dfm.base.api.timeseries;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.api.timeseries.TsPeriod;

/**
 *
 * @author Jean Palate
 */
public class TsInformationUpdates {

    /**
     *
     */
    @lombok.AllArgsConstructor
    public static class Update {

        /**
         *
         */
        @lombok.Getter
        private final TsPeriod period;
        /**
         *
         */
        @lombok.Getter
        private final int series;

        @lombok.Getter
        private final double observation;

        @lombok.Getter
        @lombok.Setter
        private double forecast;

        public Update(final TsPeriod period, final int series, final double y) {
            this.period = period;
            this.series = series;
            this.observation = y;
            this.forecast = y;
        }

        /**
         *
         * @return
         */
        public double getNews() {
            return observation - forecast;
        }

    }

    private final List<Update> news = new ArrayList<>();
    private final List<Update> revisions = new ArrayList<>();

    /**
     *
     * @param update
     */
    public void addNews(Update update) {
        news.add(update);
    }

    public void addRevision(Update update) {
        revisions.add(update);
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
                first = TsUtility.endPeriod(update.period, freq);
            } else {
                TsPeriod cur = TsUtility.endPeriod(update.period, freq);
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
                last = TsUtility.endPeriod(update.period, freq);
            } else {
                TsPeriod cur = TsUtility.endPeriod(update.period, freq);
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
    public static TsDomain updatesDomain(int freq, List<Update> updates) {
        TsPeriod first = null;
        TsPeriod last = null;
        for (Update update : updates) {
            TsPeriod cur = TsUtility.lastPeriod(update.period, freq);
            if (first == null || cur.isBefore(first)) {
                first = cur;
            }
            cur = cur.plus(1);
            if (last == null || cur.isAfter(last)) {
                last = cur;
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
    
    public List<Update> all(){
        ArrayList<Update> all=new ArrayList<>(news);
        all.addAll(revisions);
        return Collections.unmodifiableList(all);
    }
    
    public boolean isEmpty(){
        return news.isEmpty() && revisions.isEmpty();
    }
}
