/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.dfm.base.core.extractors;

import jdplus.dfm.base.api.DfmDictionaries;
import jdplus.dfm.base.core.DfmResultsNews;
import jdplus.toolkit.base.api.information.InformationExtractor;
import jdplus.toolkit.base.api.information.InformationMapping;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import nbbrd.service.ServiceProvider;

/**
 *
 * @author LEMASSO
 */
@ServiceProvider(InformationExtractor.class)
public class DfmResultsNewsExtractor extends InformationMapping<DfmResultsNews> {
    
    public DfmResultsNewsExtractor() {
        set(DfmDictionaries.SERIES_INDEX, int[].class, source -> source.getSeriesIndex());
        set(DfmDictionaries.SERIES_PERIOD, String[].class, source -> source.getSeriesPeriod());
        set(DfmDictionaries.SERIES_EXPECTED_VALUE_T, double[].class, source -> source.getSeriesExpectedValueT().toArray());
        set(DfmDictionaries.SERIES_EXPECTED_VALUE, double[].class, source -> source.getSeriesExpectedValue().toArray());
        set(DfmDictionaries.SERIES_OBSERVED_VALUE_T, double[].class, source -> source.getSeriesObservedValueT().toArray());
        set(DfmDictionaries.SERIES_OBSERVED_VALUE, double[].class, source -> source.getSeriesObservedValue().toArray());
        set(DfmDictionaries.SERIES_NEWS_T, double[].class, source -> source.getSeriesNewsT().toArray());
        set(DfmDictionaries.SERIES_NEWS, double[].class, source -> source.getSeriesNews().toArray());
        set(DfmDictionaries.SERIES_WEIGHTS_T, Matrix.class, source -> source.getSeriesWeightsT());
        set(DfmDictionaries.SERIES_WEIGHTS, Matrix.class, source -> source.getSeriesWeights());
        set(DfmDictionaries.SERIES_IMPACTS_T, Matrix.class, source -> source.getSeriesImpactsT());
        set(DfmDictionaries.SERIES_IMPACTS, Matrix.class, source -> source.getSeriesImpacts());
        set(DfmDictionaries.FORECASTS_PERIODS, String[].class, source -> source.getForecastsPeriods());
        set(DfmDictionaries.OLD_FORECASTS_T, double[].class, source -> source.getOldForecastsT().toArray());
        set(DfmDictionaries.OLD_FORECASTS, double[].class, source -> source.getOldForecasts().toArray());
        set(DfmDictionaries.REVISED_FORECASTS_T, double[].class, source -> source.getRevisedForecastsT().toArray());
        set(DfmDictionaries.REVISED_FORECASTS, double[].class, source -> source.getRevisedForecasts().toArray());
        set(DfmDictionaries.NEW_FORECASTS_T, double[].class, source -> source.getNewForecastsT().toArray());
        set(DfmDictionaries.NEW_FORECASTS, double[].class, source -> source.getNewForecasts().toArray());
    }
    
    @Override
    public Class<DfmResultsNews> getSourceClass() {
        return DfmResultsNews.class;
    }
}
