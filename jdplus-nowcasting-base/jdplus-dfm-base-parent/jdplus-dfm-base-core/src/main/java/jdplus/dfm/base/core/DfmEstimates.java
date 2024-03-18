/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.dfm.base.core;

import jdplus.toolkit.base.api.information.GenericExplorable;
import jdplus.toolkit.base.api.math.matrices.Matrix;

/**
 *
 * @author LEMASSO
 */
@lombok.Value
@lombok.Builder(builderClassName="Builder")
public class DfmEstimates implements GenericExplorable {

  private DynamicFactorModel dfm;
  private double ll;
  private Matrix hessian;
  private double[] gradient;
  private boolean hasConverged;
  
}
