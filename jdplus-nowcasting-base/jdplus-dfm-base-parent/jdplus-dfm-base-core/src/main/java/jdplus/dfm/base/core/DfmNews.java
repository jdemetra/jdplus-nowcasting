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
package jdplus.dfm.base.core;

import java.util.List;
import jdplus.dfm.base.api.timeseries.TsInformationSet;
import jdplus.dfm.base.api.timeseries.TsInformationUpdates;
import jdplus.dfm.base.api.timeseries.TsInformationUpdates.Update;
import jdplus.dfm.base.api.timeseries.TsUtility;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.timeseries.TsDomain;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.data.DataBlockIterator;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.LowerTriangularMatrix;
import jdplus.toolkit.base.core.math.matrices.SymmetricMatrix;
import jdplus.toolkit.base.core.ssf.State;
import jdplus.toolkit.base.core.ssf.StateInfo;
import jdplus.toolkit.base.core.ssf.StateStorage;
import jdplus.toolkit.base.core.ssf.multivariate.IMultivariateSsf;
import jdplus.toolkit.base.core.ssf.multivariate.MultivariateFilteringInformation;
import jdplus.toolkit.base.core.ssf.multivariate.MultivariateOrdinaryFilter;
import jdplus.toolkit.base.core.ssf.multivariate.MultivariateOrdinarySmoother;
import jdplus.toolkit.base.core.ssf.multivariate.SsfMatrix;

/**
 * Computation of the news (See Banbura an Modagno, appendix D of the reference
 * paper for further details) This implementation considers separately new
 * figures and revised figures, to get "manageable" state space models
 *
 * @author Jean Palate
 */
public class DfmNews {

    private final DynamicFactorModel model;
    private final IMultivariateSsf ssf;

    /**
     * states relative to the old data (excluded missing values). Computed on
     * [oldDomain.end, fullDomain.end[ (using the indices of the full domain)
     */
    private StateStorage oldStates;
    /**
     * states relative to the new data (excluded missing values). Computed on
     * [newDomain.end, fullDomain.end[ (using the indices of the full domain)
     */
    private StateStorage newStates;
    /**
     * states relative to the revised data (excluded missing values). Computed
     * on [oldDomain.start, fullDomain.end[ (using the indices of the full
     * domain)
     */
    private StateStorage revisedStates;

    /**
     * Number of lags used to compute the impact of the news/revisions
     */
    private int nbnews, nbrev;

    private TsInformationSet oldSet, newSet, revisedSet;

    private TsInformationUpdates updates;
    private FastMatrix covNews, lcovNews;
    private FastMatrix covRevisions, lcovRevisions;

    private int ext = 2;
    /**
     * common domain (= domain where all the series are defined) of the old/new
     * data set
     */
    private TsDomain oldDomain, newDomain;
    /**
     * "Complete" domain, which is the union of the domains of the old
     * information set and of the new information set
     */
    private TsDomain fullDomain;
    /**
     * domain containing all the news
     */
    private TsDomain newsDomain;
    /**
     * domain containing all the revisions
     */
    private TsDomain revisionsDomain;

    /**
     *
     * @param model
     */
    public DfmNews(DynamicFactorModel model) {
        this.model = model;
        this.ssf = model.ssfRepresentation(0);
    }

    public TsInformationSet getOldInformationSet() {
        return this.oldSet;
    }

    public TsInformationSet getNewInformationSet() {
        return this.newSet;
    }

    /**
     * Returns the old data set with revisions from the new data set. The news
     * are not present.
     *
     * @return Old data + revisions
     */
    public TsInformationSet getRevisedInformationSet() {
        return this.revisedSet;
    }

    public DynamicFactorModel getModel() {
        return model;
    }

    public TsDomain getNewsDomain() {
        return newsDomain;
    }

    public TsDomain getRevisionsDomain() {
        return revisionsDomain;
    }

    /**
     * Computes the news between two consecutive information sets
     *
     * @param oldSet The old information set
     * @param newSet The new information set
     * @return True if the news have been successfully computed
     */
    public boolean process(TsInformationSet oldSet, TsInformationSet newSet) {
        this.oldSet = oldSet;
        this.newSet = newSet;
        // recall : revisions = revisedSet - oldSet, news = newSet - revisedSet
        revisedSet = this.oldSet.revisedData(this.newSet);
        updates = this.oldSet.updates(this.newSet);    // Calculates news and revisions
        if (updates.isEmpty()) {
            return false;
        }
        computeDomains();
        return calcNews();
    }

    private boolean calcNews() {
        // Calculates News
        FastMatrix M_old = FastMatrix.of(oldSet.generateMatrix(fullDomain));
        FastMatrix M_rev = FastMatrix.of(revisedSet.generateMatrix(fullDomain));
        FastMatrix M_new = FastMatrix.of(newSet.generateMatrix(fullDomain));

        oldStates = smoothData(M_old);
        if (oldStates == null) {
            return false;
        }

        // Calculates Revisions
        if (!updates.revisions().isEmpty()) {
            revisedStates = smoothData(M_rev);
            if (revisedStates == null) {
                return false;
            }
//            computeRevisionsCovariance(M_old);
        } else {
            revisedStates = oldStates;
        }

        updateNews();

        if (!updates.news().isEmpty()) {
            newStates = smoothData(M_new);
            if (newStates == null) {
                return false;
            }
            computeNewsCovariance(M_rev);
        } else {
            newStates = revisedStates;
        }

        return true;
    }

    public double getOldForecast(int series, TsPeriod p) {
        int pos = fullDomain.getStartPeriod().until(TsUtility.lastPeriod(p, fullDomain.getAnnualFrequency()));
        DataBlock A = oldStates.a(pos);
        return A != null ? ssf.loading(series).ZX(pos, A) : Double.NaN;
    }

    public double getRevisedForecast(int series, TsPeriod p) {
        int pos = fullDomain.getStartPeriod().until(TsUtility.lastPeriod(p, fullDomain.getAnnualFrequency()));
        DataBlock A = revisedStates.a(pos);
        return A != null ? ssf.loading(series).ZX(pos, A) : Double.NaN;
    }

    public double getNewForecast(int series, TsPeriod p) {
        int pos = fullDomain.getStartPeriod().until(TsUtility.lastPeriod(p, fullDomain.getAnnualFrequency()));
        DataBlock A = newStates.a(pos);
        return A != null ? ssf.loading(series).ZX(pos, A) : Double.NaN;
    }

    public int getMaxNewsExtensionPeriod() {
        return ext;
    }

    public void setMaxNewsExtensionPeriod(int n) {
        if (n != ext) {
            ext = n;
            calcNews();
        }
    }

    private void computeDomains() {
        TsDomain d0 = revisedSet.getCurrentDomain(); // same as oldSet
        TsDomain d1 = newSet.getCurrentDomain();
        fullDomain = d0.union(d1);
        // Intersection of the domains, after excluding missings at the extremities
        oldDomain = revisedSet.actualData().getCommonDomain();
        newDomain = newSet.actualData().getCommonDomain();
        computeNewsDomain();
        computeRevisionsDomain();
    }

    /**
     * The news domain is defined as the domain starting at the end of the
     * common domain of the old data. So, it starts at the first missing value
     * (excluding starting/intermediary missing values) and it ends at the last
     * news
     */
    private void computeNewsDomain() {
        int freq = fullDomain.getAnnualFrequency();
        if (updates.news().isEmpty()) {
            newsDomain = null;
            nbnews = 0;
        } else {
            newsDomain = TsInformationUpdates.updatesDomain(freq, updates.news());
//            TsDomain ndomain = TsInformationUpdates.updatesDomain(freq, updates.news());
//            TsPeriod start = oldDomain.getEndPeriod();
//            newsDomain = TsDomain.of(start, start.until(ndomain.getEndPeriod()));
            int nb = model.defaultSsfBlockLength();
            // number of items we have to include in the state space for each factor 
            // if t0 = newsDomain.start, t1 = fullDomain.last and m = measurementLength
            // we need to have items in [t0-m+1, t1[ (and more that the default block length) 
            nbnews = Math.max(nb, model.measurementsLength() + newsDomain.getStartPeriod().until(fullDomain.getEndPeriod()));
        }
    }

    /**
     * For the revisions domain, we go back to 'ext' years of revisions
     */
    private void computeRevisionsDomain() {
        int freq = fullDomain.getAnnualFrequency();
        if (updates.revisions().isEmpty()) {
            revisionsDomain = null;
            nbrev = 0;
        } else {
            TsDomain rdomain = TsInformationUpdates.updatesDomain(freq, updates.revisions());
            TsPeriod start = rdomain.getStartPeriod();
            int n = start.until(oldDomain.getEndPeriod());
            if (n > ext * freq) {
                rdomain = rdomain.drop(n - ext * freq, 0);
            }
            revisionsDomain = rdomain;
            int nb = model.defaultSsfBlockLength();
            // number of items we have to include in the state space for each factor 
            // if t0 = newsDomain.start, t1 = fullDomain.last and m = measurementLength
            // we need to have items in [t0-m+1, t1[ (and more that the default block length) 
            nbrev = Math.max(nb, model.measurementsLength() + revisionsDomain.getStartPeriod().until(fullDomain.getEndPeriod()));
        }
    }

    /**
     *
     * @return
     */
    private StateStorage smoothData(FastMatrix M) {
        // We don't compute the variances
        MultivariateOrdinarySmoother smoother = MultivariateOrdinarySmoother.builder(ssf)
                .calcVariance(false)
                .calcSmoothationsVariance(false)
                .build();
        SsfMatrix data = new SsfMatrix(M);
        StateStorage ss = StateStorage.light(StateInfo.Smoothed);
        MultivariateOrdinaryFilter filter = new MultivariateOrdinaryFilter();
        MultivariateFilteringInformation fresults = new MultivariateFilteringInformation();
        if (filter.process(ssf, data, fresults)) {
            ss.prepare(ssf.getStateDim(), 0, data.getObsCount());
            smoother.process(0, data.getObsCount(), fresults, ss);
            return ss;
        } else {
            return null;
        }
    }

    private StateStorage smoothDataEx(int nb, FastMatrix M, TsPeriod start) {
        IMultivariateSsf xssf = model.ssfRepresentationWithBlockLength(nb);
        MultivariateOrdinarySmoother smoother = MultivariateOrdinarySmoother.builder(xssf)
                .calcVariance(true)
                .calcSmoothationsVariance(true)
                .build();
        int n = fullDomain.getStartPeriod().until(start);
        SsfMatrix data = new SsfMatrix(M);
        StateStorage ss = StateStorage.full(StateInfo.Smoothed);
        MultivariateOrdinaryFilter filter = new MultivariateOrdinaryFilter();
        MultivariateFilteringInformation fresults = new MultivariateFilteringInformation();
        if (filter.process(xssf, data, fresults)) {
            ss.prepare(xssf.getStateDim(), n, data.getObsCount());
            smoother.process(n, data.getObsCount(), fresults, ss);
            return ss;
        } else {
            return null;
        }
    }

    /**
     * Updates the news with the forecasts computed on the old data
     */
    private void updateNews() {
        int freq = fullDomain.getAnnualFrequency();
        TsPeriod start = fullDomain.getStartPeriod();
        for (Update update : updates.news()) {
            int pos = start.until(TsUtility.lastPeriod(update.getPeriod(), freq));
            update.setForecast(ssf.loading(update.getSeries()).ZX(pos, revisedStates.a(pos)));
        }
    }

    /**
     *
     * @return
     */
    public StateStorage getOldSmoothingResults() {
        return oldStates;
    }

    /**
     *
     * @return
     */
    public StateStorage getNewSmoothingResults() {
        return newStates;
    }

    private void computeNewsCovariance(FastMatrix M) {

        StateStorage ss = smoothDataEx(nbnews, M, newsDomain.getStartPeriod());
        covNews = ss.P(fullDomain.length() - 1).deepClone();

        int freq = fullDomain.getAnnualFrequency();
        List<Update> lupdates = this.updates.news();
        int nupdates = lupdates.size();
        int c = model.measurementsLength();
        int nf = model.getNfactors();
        int d = c * nf;
        int n = fullDomain.length() - 1;
        lcovNews = FastMatrix.square(nupdates);
        FastMatrix V = FastMatrix.square(d);
        DataBlockIterator vcols = V.columnsIterator();
        DataBlock tmp = DataBlock.make(d);
        TsPeriod end = fullDomain.getEndPeriod();
        for (int i = 0; i < nupdates; ++i) {
            Update iupdate = lupdates.get(i);
            int istart = TsUtility.endPeriod(iupdate.getPeriod(), freq).until(end);
            for (int j = 0; j <= i; ++j) {
                Update jupdate = lupdates.get(j);
                // copy the right covariance
                int jstart = TsUtility.endPeriod(jupdate.getPeriod(), freq).until(end);
                V.set(0);
                for (int r = 0; r < nf; ++r) {
                    for (int s = 0; s < nf; ++s) {
                        V.extract(r * c, c, s * c, c).copy(covNews.extract(r * nbnews + istart, c,
                                s * nbnews + jstart, c));
                    }
                }

                vcols.begin();
                tmp.set(vcols, col -> ssf.loading(iupdate.getSeries()).ZX(n - istart, col));
                double q = ssf.loading(jupdate.getSeries()).ZX(n - jstart, tmp);
                if (i == j) {
                    q += model.getMeasurements().get(iupdate.getSeries()).getVariance();
                }
                lcovNews.set(i, j, q);
            }
        }
        SymmetricMatrix.fromLower(lcovNews);
        SymmetricMatrix.lcholesky(lcovNews, State.ZERO);
        LowerTriangularMatrix.toLower(lcovNews);
    }

//    private void computeRevisionsCovariance(FastMatrix M) {
//        StateStorage ss = smoothDataEx(nbrev, M, revisionsDomain.getStartPeriod());
//        covRevisions = ss.P(fullDomain.length() - 1).deepClone();
//
//        int freq = fullDomain.getAnnualFrequency();
//        List<Update> lupdates = updates.revisions();
//        int nupdates = lupdates.size();
//        int c = model.measurementsLength();
//        int nf = model.getNfactors();
//        int d = c * nf;
//        int n = fullDomain.length() - 1;
//        lcovRevisions = FastMatrix.square(nupdates);
//        FastMatrix V = FastMatrix.square(d);
//        DataBlockIterator vcols = V.columnsIterator();
//        DataBlock tmp = DataBlock.make(d);
//        TsPeriod end = fullDomain.getEndPeriod();
//        for (int i = 0; i < nupdates; ++i) {
//            Update iupdate = lupdates.get(i);
//            int istart = TsUtility.endPeriod(iupdate.getPeriod(), freq).until(end);
//            for (int j = 0; j <= i; ++j) {
//                Update jupdate = lupdates.get(j);
//                // copy the right covariance
//                int jstart = TsUtility.endPeriod(jupdate.getPeriod(), freq).until(end);
//                V.set(0);
//                for (int r = 0; r < nf; ++r) {
//                    for (int s = 0; s < nf; ++s) {
//                        V.extract(r * c, c, s * c, c).copy(covRevisions.extract(r * nbrev + istart, c,
//                                s * nbrev + jstart, c));
//                    }
//                }
//
//                vcols.begin();
//                tmp.set(vcols, col -> ssf.loading(iupdate.getSeries()).ZX(n - istart, col));
//                double q = ssf.loading(jupdate.getSeries()).ZX(n - jstart, tmp);
//                if (i == j) {
//                    q += model.getMeasurements().get(iupdate.getSeries()).getVariance();
//                }
//                lcovRevisions.set(i, j, q);
//            }
//        }
//        SymmetricMatrix.fromLower(lcovRevisions);
//        SymmetricMatrix.lcholesky(lcovRevisions, State.ZERO);
//        LowerTriangularMatrix.toLower(lcovRevisions);
//    }
//
    /**
     *
     * @return
     */
    public FastMatrix getStateCovariance() {
        return covNews;
    }

    /**
     *
     * @return
     */
    public TsInformationUpdates newsDetails() {
        return updates;
    }

    /**
     *
     * @return
     */
    public DoubleSeq news() {
        List<Update> news = updates.news();
        int n = news.size();
        double[] a = new double[n];
        int i = 0;
        for (Update cnews : news) {
            a[i++] = cnews.getNews();
        }
        return DoubleSeq.of(a);
    }

    public DoubleSeq revisions() {
        List<Update> revisions = updates.revisions();
        int n = revisions.size();
        double[] a = new double[n];
        int i = 0;
        for (Update rev : revisions) {
            a[i++] = rev.getNews();
        }
        return DoubleSeq.of(a);
    }

    public DoubleSeq weights(int series, TsPeriod p) {
        List<Update> lupdates = this.updates.news();
        if (lupdates.isEmpty()) {
            return DoubleSeq.empty();
        }
        int nupdates = lupdates.size();
        DataBlock a = DataBlock.make(nupdates);
        int nf = model.getNfactors();
        int c = model.measurementsLength();
        int d = c * nf;
        int n = fullDomain.length() - 1;
        FastMatrix V = FastMatrix.square(d);
        DataBlockIterator vcols = V.columnsIterator();
        DataBlock tmp = DataBlock.make(d);
        TsPeriod end = fullDomain.getEndPeriod();
        int freq = end.annualFrequency();
        int istart = TsUtility.endPeriod(p, freq).until(end);
        for (int j = 0; j < nupdates; ++j) {
            Update jupdate = lupdates.get(j);
            int jstart = TsUtility.endPeriod(jupdate.getPeriod(), freq).until(end);
            V.set(0);
            for (int r = 0; r < nf; ++r) {
                for (int s = 0; s < nf; ++s) {
                    V.extract(r * c, c, s * c, c).copy(covNews.extract(r * nbnews + istart, c,
                            s * nbnews + jstart, c));
                }
            }
            vcols.begin();
            tmp.set(vcols, col -> ssf.loading(series).ZX(n - istart, col));
            double q = ssf.loading(jupdate.getSeries()).ZX(n - jstart, tmp);
            a.set(j, q);
        }
        // w = A * (LL')^-1 <-> w(LL')=A
        // B = wL, BL' = A <-> LB'=A'
        LowerTriangularMatrix.solveLx(lcovNews, a, State.ZERO); // B
        LowerTriangularMatrix.solvexL(lcovNews, a, State.ZERO);
        return a;
    }

//    public DoubleSeq weightsRevisions(int series, TsPeriod p) {
//        List<Update> lupdates = this.updates.revisions();
//        if (lupdates.isEmpty()) {
//            return DoubleSeq.empty();
//        }
//        int nupdates = lupdates.size();
//        DataBlock a = DataBlock.make(nupdates);
//        int c = model.measurementsLength();
//        int nf = model.getNfactors();
//        int d = c * nf;
//        int n = fullDomain.length() - 1;
//        FastMatrix V = FastMatrix.square(d);
//        DataBlockIterator vcols = V.columnsIterator();
//        DataBlock tmp = DataBlock.make(d);
//        TsPeriod end = fullDomain.getEndPeriod();
//        int freq = end.annualFrequency();
//        int istart = TsUtility.endPeriod(p, freq).until(end);
//        for (int j = 0; j < nupdates; ++j) {
//            Update jupdate = lupdates.get(j);
//            int jstart = TsUtility.endPeriod(jupdate.getPeriod(), freq).until(end);
//            V.set(0);
//            for (int r = 0; r < nf; ++r) {
//                for (int s = 0; s < nf; ++s) {
//                    V.extract(r * c, c, s * c, c).copy(covRevisions.extract(r * nbrev + istart, c,
//                            s * nbrev + jstart, c));
//                }
//            }
//            tmp.set(0);
//            vcols.begin();
//            tmp.set(vcols, col -> ssf.loading(series).ZX(n - istart, col));
//            double q = ssf.loading(jupdate.getSeries()).ZX(n - jstart, tmp);
//            a.set(j, q);
//        }
//        // w = A * (LL')^-1 <-> w(LL')=A
//        // B = wL, BL' = A <-> LB'=A'
//        LowerTriangularMatrix.solveLx(lcovRevisions, a, State.ZERO); // B
//        LowerTriangularMatrix.solvexL(lcovRevisions, a, State.ZERO);
//        return a;
//    }
}
