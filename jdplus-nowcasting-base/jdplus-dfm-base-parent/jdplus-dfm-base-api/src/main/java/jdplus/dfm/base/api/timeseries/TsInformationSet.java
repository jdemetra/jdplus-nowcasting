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

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import jdplus.dfm.base.api.timeseries.TsInformationUpdates.Update;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.timeseries.TimeSelector;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsDataTable;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.api.timeseries.TsUnit;
import jdplus.toolkit.base.api.util.IntList;

/**
 *
 * @author Jean Palate
 */
public class TsInformationSet {

    /**
     *
     * @param input
     */
    public TsInformationSet(Collection<TsData> input) {
        table = TsDataTable.of(input);
    }

    public TsData[] toArray() {
        return table.getData().toArray(TsData[]::new);
    }

    /**
     * Creates a new information set with only the revised data in comparison
     * with this data set (the domains of the series of the returned data set
     * are identical to the domains of the series of this information set.
     * Missing obs in this data set are also considered as missing in the
     * revised data (so, we consider actual revisions)
     *
     * So, the revisedData contains the new figures with the same "missing
     * pattern" as the old one.
     *
     * @param newdata
     * @return
     */
    public TsInformationSet revisedData(TsInformationSet newdata) {
        List<TsData> ndata = new ArrayList<>();
        int n = getSeriesCount();
        for (int i = 0; i < n; ++i) {
            TsData cur = series(i);
            TsData ncur = newdata.series(i);
            ncur = TsData.fitToDomain(ncur, cur.getDomain());
            double[] nvals = null;
            for (int j = 0; j < cur.length(); ++j) {
                if (TsUtility.isMissing(cur, j) && !TsUtility.isMissing(ncur, j)) {
                    if (nvals == null) {
                        nvals = ncur.getValues().toArray();
                    }
                    nvals[j] = Double.NaN;
                }
            }
            if (nvals != null) {
                ncur = TsData.ofInternal(ncur.getStart(), nvals);
            }
            ndata.add(ncur);
        }
        return new TsInformationSet(ndata);
    }

    public TsInformationSet actualData() {
        return new TsInformationSet(table.getData()
                .stream()
                .peek(s -> s.cleanExtremities())
                .toList());
    }

    public List<TsData> all() {
        return table.getData();
    }

    /**
     * Extends the series so that the last periods of the series are the last
     * complete period strictly before the given date
     *
     * @param end End (excluded)
     * @return
     */
    public TsInformationSet extendTo(final LocalDate end) {
        List<TsData> ndata = new ArrayList<>();
        for (TsData s : table.getData()) {
            ndata.add(TsUtility.extend(s, end));
        }
        return new TsInformationSet(ndata);
    }

    /**
     * Union of all the domains (at the highest periodicity)
     *
     * @return
     */
    public TsDomain getCurrentDomain() {
        return table.getDomain();
    }

    /**
     * Intersection of all the domains (at the highest periodicity)
     *
     * @return
     */
    public TsDomain getCommonDomain() {
        if (table.getData().isEmpty()) {
            return TsDomain.DEFAULT_EMPTY;
        }
        int f = table.getDomain().getAnnualFrequency();
        TsDomain common = null;
        int n = getSeriesCount();
        for (int i = 0; i < n; ++i) {
            TsDomain cur = series(i).getDomain();
            TsPeriod p0 = TsPeriod.of(TsUnit.ofAnnualFrequency(f), cur.start());
            TsPeriod p1 = TsPeriod.of(TsUnit.ofAnnualFrequency(f), cur.end());
            TsDomain fcur = TsDomain.of(p0, p0.until(p1));
            common = common != null ? common.intersection(fcur) : fcur;
        }
        return common;
    }

    /**
     * Number of series
     *
     * @return
     */
    public int getSeriesCount() {
        return table.getData().size();
    }

    /**
     * Number of actual figures (non missing)
     *
     * @return
     */
    public int getDataCount() {
        return table.getData()
                .stream()
                .mapToInt(s -> s.getValues().count(x -> Double.isFinite(x)))
                .sum();
    }

    /**
     * Retrieve a specific series
     *
     * @param idx
     * @return
     */
    public TsData series(int idx) {
        return table.getData().get(idx);
    }

    /**
     * Gets the data corresponding to the given domain or to the union domain if
     * the given domain is null
     *
     * @param domain
     * @return
     */
    public Matrix generateMatrix(final TsDomain domain) {
        TsDomain tdomain = table.getDomain();
        if (tdomain.isEmpty()) {
            return null;
        }
        if (domain == null) {
            return generateMatrix(tdomain);
        }
        if (domain.getAnnualFrequency() != tdomain.getAnnualFrequency()) {
            return null;
        }

        int nr = domain.length(), nc = table.getData().size();
        double[] data = new double[nr * nc];
        Arrays.fill(data, Double.NaN);

        TsDomain common = tdomain.intersection(domain);
        int r0 = domain.getStartPeriod().until(common.getStartPeriod()); // del relative to the requested domain
        int tr0 = tdomain.getStartPeriod().until(common.getStartPeriod()); // del relative to th current domain
        int trmax = tr0 + common.length();
        TsDataTable.Cursor cursor = table.cursor(TsDataTable.DistributionType.LAST);
        for (int i = tr0, j = r0; i < trmax; ++i, ++j) {
            for (int c = 0; c < nc; ++c) {
                cursor.moveTo(i, c);
                TsDataTable.ValueStatus status = cursor.getStatus();
                if (status == TsDataTable.ValueStatus.PRESENT) {
                    data[j + c * nr] = cursor.getValue();
                }
            }
        }
        return Matrix.of(data, nr, nc);
    }

    /**
     * Fill in periods for each series where new data are present (does not take
     * into account values that have been revised). Takes new data before first
     * element in old dataset, new data after last element in old dataset and
     * new data where there was a missing value and data revisions. Both types
     * of updates are put in separate lists. At this stage, the updates don't
     * contain forecasts
     *
     * @param ndata New data
     * @return List of newly added data
     */
    public TsInformationUpdates updates(TsInformationSet ndata) {
        int n = getSeriesCount();
        if (n != ndata.getSeriesCount()) {
            return null;
        }
        TsInformationUpdates updates = new TsInformationUpdates();
        for (int i = 0; i < n; ++i) {
            TsData olds = series(i), news = ndata.series(i);
            TsPeriod start = news.getStart();
            int del = olds.getStart().until(start);
            for (int j = 0; j < news.length(); ++j) {
                double y = news.getValue(j);
                if (!Double.isNaN(y)) {
                    int k = j + del;
                    if (k < 0 || k >= olds.length() || TsUtility.isMissing(olds, k)) {
                        updates.addNews(new Update(start.plus(j), i, y));
                    }
                }
            }

            // Calculates revisions
            start = olds.getStart();
            TsData newFit = TsData.fitToDomain(news, olds.getDomain());
            for (int j = 0; j < olds.length(); ++j) {
                double y = newFit.getValue(j), fy = olds.getValue(j);
                if (!Double.isNaN(y)
                        && !Double.isNaN(fy)
                        && y != olds.getValue(j)) {
                    updates.addRevision(new Update(start.plus(j), i, y, fy));
                }
            }
        }
        return updates;
    }

    public LocalDate[] generatePublicationCalendar(int[] delays) {
        SortedSet<LocalDate> sdays = new TreeSet<>();
        List<TsData> allSeries = table.getData();
        int i = 0;
        for (TsData s : allSeries) {
            TsDomain dom = s.getDomain();
            int ndel = delays == null ? 0 : delays[i++];
            for (int j = 0; j < s.length(); ++j) {
                if (!TsUtility.isMissing(s, j)) {
                    TsPeriod p = dom.get(j);
                    // ndel days after the end of the series. 1 day more than in the first implementation !
                    // if ndel is 0, the publication day is the first day after the observed periods
                    LocalDate pub = p.end().toLocalDate().plusDays(ndel);
                    sdays.add(pub);
                }
            }
        }
        LocalDate[] days = new LocalDate[sdays.size()];
        return sdays.toArray(days);
    }

    public LocalDate[] generatePublicationCalendar(IntList delays, LocalDate start) {
        SortedSet<LocalDate> sdays = new TreeSet<>();
        int n = getSeriesCount();
        for (int i = 0; i < n; ++i) {
            TsData s = series(i);
            TsDomain dom = s.getDomain();
            int ndel = (delays == null || delays.isEmpty()) ? 0 : delays.get(i);
            int pos = dom.indexOf(start.atStartOfDay());
            if (pos < 0 && start.isBefore(dom.start().toLocalDate())) {
                pos = 0;
            }
            if (pos >= 0) {
                for (int j = pos; j < s.length(); ++j) {
                    if (!TsUtility.isMissing(s, j)) {
                        TsPeriod p = dom.get(j);
                        LocalDate pub = p.end().toLocalDate().plusDays(ndel);
                        sdays.add(pub);
                    }
                }
            }
        }
        LocalDate[] days = new LocalDate[sdays.size()];
        return sdays.toArray(days);
    }

    /**
     *
     * @param delays
     * @param date
     * @return
     */
    public TsInformationSet generateInformation(final IntList delays, final LocalDate date) {
        int n = getSeriesCount();
        List<TsData> inputc = new ArrayList<>(n);
        for (int i = 0; i < n; ++i) {
            LocalDate last = date;
            if (delays != null && !delays.isEmpty() && i < delays.size()) {
                last = last.minusDays(delays.get(i));
            }
            TimeSelector sel = TimeSelector.to(last.atStartOfDay());
            inputc.add(series(i).select(sel));
        }
        return new TsInformationSet(inputc);
    }

    private final TsDataTable table;
}
