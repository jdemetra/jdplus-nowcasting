/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.dfm.base.core;

import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.information.GenericExplorable;
import jdplus.toolkit.base.api.math.matrices.Matrix;

/**
 *
 * @author LEMASSO
 */

@lombok.Value
@lombok.Builder(builderClassName = "Builder")
public class DfmResultsNews implements GenericExplorable{
    
    int[] seriesIndex;
    String[] seriesPeriod;
    DoubleSeq seriesExpectedValueT; 
    DoubleSeq seriesExpectedValue;   
    DoubleSeq seriesObservedValueT;
    DoubleSeq seriesObservedValue;
    DoubleSeq seriesNewsT;
    DoubleSeq seriesNews;
    Matrix seriesWeightsT;
    Matrix seriesWeights;
    Matrix seriesImpactsT;
    Matrix seriesImpacts;
    String[] forecastsPeriods;
    DoubleSeq oldForecastsT;
    DoubleSeq oldForecasts;
    DoubleSeq revisedForecastsT; // if no revision, revised = old
    DoubleSeq revisedForecasts;
    DoubleSeq newForecastsT;
    DoubleSeq newForecasts;
}
