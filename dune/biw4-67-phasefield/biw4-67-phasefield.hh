#ifndef BIW4_67_PHASEFIELD_HH
#define BIW4_67_PHASEFIELD_HH


#include <dune/common/fvector.hh> 
#include <dune/common/fmatrix.hh>
// add your classes here


namespace Dune {
    namespace Biw467 {
        

        template<int dim, class field_type=double>
        class material {

            public:

            double _griffithReleaseRate,_shearModulus,_firstLameConstant,_regularisationParameter;

            material(double griffithReleaseRate,double shearModulus,double firstLameConstant,double regularisationParameter) :
                _griffithReleaseRate(griffithReleaseRate),
                _shearModulus(shearModulus),
                _firstLameConstant(firstLameConstant),
                _regularisationParameter(regularisationParameter) {};

            double phasefieldEnergyDensity(double phasefieldNumber, const Dune::FieldVector<field_type, dim> gradient) const {
                double normgrad = gradient.two_norm2();
                return _griffithReleaseRate * 1.0/(2.0 * _regularisationParameter)*(phasefieldNumber*phasefieldNumber + _regularisationParameter*_regularisationParameter*normgrad);
            }

            double degradationFunction(double phasefieldNumber) const {
                return (1.0 - phasefieldNumber)*(1.0 - phasefieldNumber);
            }

            double degradationFunctionDerivative(double phasefieldNumber) const {
                return -2.0*(1.0 - phasefieldNumber);
            }

            double strainEnergyDensity(const Dune::FieldMatrix<field_type, dim, dim>& strains) const {
                
                double trace = 0.0;
                double traceE2 = 0.0;
                for (int i = 0; i <dim; i++) {
                    trace += strains[i][i];
                    for (int j = 0; j < dim; j++) {
                        traceE2 += strains[i][j] * strains[j][i];
                    }
                }
                

                return _shearModulus*traceE2 + 0.5*_firstLameConstant*trace*trace;
            }

            void stresses(const Dune::FieldMatrix<field_type, dim, dim>& strains, Dune::FieldMatrix<field_type, dim, dim>& stresses) const {
                stresses = 0.0;
                double trace = 0.0;
                for (int i = 0; i <dim; i++) {
                    trace += strains[i][i];
                }
                for (int i = 0; i < dim; i++) {
                    stresses[i][i] +=  _firstLameConstant * trace;
                    for (int j; j < dim; j++) {
                        stresses[i][j] += 2.0 * _shearModulus * strains[i][j];
                    }
                }
            }

        };
    }
}


#endif // BIW4_67_PHASEFIELD_HH
