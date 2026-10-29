# Copyright (C) 2026 Zachary Kenneth Stewart

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import math
import numpy as np

class BinaryTrait:
    '''
    Handling of a binary trait: its PLINK coding (1 = control, 2 = case), the
    SAIGE settings for testing it (a logistic mixed model, with SPA and Firth
    correction for P < 0.05), and its summaries in lead variant inspections
    (allele frequencies and case rates of carriers and non-carriers).

    Example:
        trait = get_trait("binary")
        value = trait.recode("2") # '1'
    '''
    TYPE = "binary"
    PHENO_NAME = "y_binary"
    MISSING_VALUES = {"-9", "0", "NA", "na", "nan", "NaN", "."}
    SAIGE_RESULT_COLUMNS = ["AF_case", "AF_ctrl"]
    LEAD_DIRECTION_COLUMNS = ["AF_CASE", "AF_CTRL", "ENRICHED_IN"]
    LEAD_CARRIER_COLUMNS = ["CARRIER_CASE_RATE", "NONCARRIER_CASE_RATE"]
    FAMILY_COLUMNS = ["N_CASES"]
    FAMILY_CARRIER_COLUMNS = ["CARRIER_CASES", "CARRIER_CASE_RATE", "NONCARRIER_CASE_RATE"]

    @property
    def settings(self):
        '''
        Trait settings recorded with a run.
        '''
        return {}

    def recode(self, value):
        '''
        Recodes a PLINK binary phenotype (1/2) to SAIGE's 0/1, or None if missing.
        '''
        if value in BinaryTrait.MISSING_VALUES:
            return None
        if value not in {"1", "2"}:
            raise ValueError(f"Phenotype value '{value}' is not a binary PLINK phenotype (1 = control, 2 = case)")
        return str(int(value) - 1)

    def summarise(self, values):
        '''
        Summarises recoded phenotype values (strings) of the samples in a model.
        '''
        cases = sum( int(v) for v in values )
        return { "cases": cases, "controls": len(values) - cases }

    @staticmethod
    def describe(summary):
        return f"{summary['cases']} cases, {summary['controls']} controls"

    def null_args(self):
        '''
        Trait arguments for SAIGE step1_fitNULLGLMM.R.
        '''
        return ["--traitType=binary"]

    def test_args(self):
        '''
        Trait arguments for SAIGE step2_SPAtests.R.
        '''
        return ["--is_Firth_beta=TRUE", "--pCutoffforFirth=0.05"]

    @staticmethod
    def rate(isCase):
        return f"{isCase.mean():.3f}" if len(isCase) > 0 else "."

    def lead_direction(self, stats, phenotypes, calls):
        '''
        Returns the LEAD_DIRECTION_COLUMNS values of a lead variant.

        Parameters:
            stats -- a dictionary of the lead's SAIGE results (see parse_saige_leads())
            phenotypes -- a numpy array of the tested samples' phenotypes
            calls -- a numpy array of their effect allele hard calls, NaN if missing
        '''
        return [stats["AF_case"], stats["AF_ctrl"], "cases" if float(stats["BETA"]) > 0 else "controls"]

    def carrier_values(self, phenotypes, isCarrier, hasCall):
        '''
        Returns the LEAD_CARRIER_COLUMNS values of a lead variant.
        '''
        isCase = phenotypes == 1
        return [BinaryTrait.rate(isCase[isCarrier]), BinaryTrait.rate(isCase[hasCall & ~isCarrier])]

    def family_values(self, phenotypes):
        '''
        Returns the FAMILY_COLUMNS values of a family's tested samples.
        '''
        return [str(int((phenotypes == 1).sum()))]

    def family_carrier_values(self, phenotypes, isCarrier, hasCall):
        '''
        Returns the FAMILY_CARRIER_COLUMNS values of a family at a lead variant.
        '''
        isCase = phenotypes == 1
        return [str(int((isCarrier & isCase).sum()))] + self.carrier_values(phenotypes, isCarrier, hasCall)

class QuantitativeTrait:
    '''
    Handling of a quantitative trait: its PLINK coding (numbers, with -9 as
    missing), the SAIGE settings for testing it (a linear mixed model, by
    default of the inverse normal transformed trait), and its summaries in lead
    variant inspections (mean phenotypes by genotype, and of carriers and
    non-carriers). Summaries are of the untransformed trait.

    Example:
        trait = get_trait("quantitative", invNormalize=True)
        value = trait.recode("3.25") # '3.25'
    '''
    TYPE = "quantitative"
    PHENO_NAME = "y_quantitative"
    MISSING_VALUES = {"-9", "NA", "na", "nan", "NaN", "."}
    SAIGE_RESULT_COLUMNS = []
    LEAD_DIRECTION_COLUMNS = ["MEAN_GENO0", "MEAN_GENO1", "MEAN_GENO2", "EFFECT"]
    LEAD_CARRIER_COLUMNS = ["CARRIER_MEAN", "NONCARRIER_MEAN"]
    FAMILY_COLUMNS = ["MEAN"]
    FAMILY_CARRIER_COLUMNS = ["CARRIER_MEAN", "NONCARRIER_MEAN"]

    def __init__(self, invNormalize=True):
        self.invNormalize = invNormalize

    @property
    def invNormalize(self):
        return self._invNormalize

    @invNormalize.setter
    def invNormalize(self, value):
        '''
        Whether SAIGE inverse normal transforms the trait before fitting the null model.
        '''
        if not isinstance(value, bool):
            raise TypeError("invNormalize must be a boolean")
        self._invNormalize = value

    @property
    def settings(self):
        return { "invNormalize": self.invNormalize }

    def recode(self, value):
        '''
        Validates a quantitative phenotype, returning it unchanged (to preserve
        its precision), or None if missing.
        '''
        if value in QuantitativeTrait.MISSING_VALUES:
            return None
        try:
            number = float(value)
        except ValueError:
            raise ValueError(f"Phenotype value '{value}' is not a number")
        if not math.isfinite(number):
            raise ValueError(f"Phenotype value '{value}' is not a finite number")
        return value

    def summarise(self, values):
        numbers = np.array([ float(v) for v in values ])
        return { "mean": round(float(numbers.mean()), 6), "sd": round(float(numbers.std(ddof=1)), 6)
                 if len(numbers) > 1 else None, "min": float(numbers.min()), "max": float(numbers.max()) }

    @staticmethod
    def describe(summary):
        sd = f"{summary['sd']:.4g}" if summary["sd"] is not None else "NA"
        return f"mean {summary['mean']:.4g}, SD {sd}, range {summary['min']:g} to {summary['max']:g}"

    def null_args(self):
        return ["--traitType=quantitative", f"--invNormalize={'TRUE' if self.invNormalize else 'FALSE'}"]

    def test_args(self):
        return []

    @staticmethod
    def mean(values):
        return f"{values.mean():.4g}" if len(values) > 0 else "."

    def lead_direction(self, stats, phenotypes, calls):
        return [ QuantitativeTrait.mean(phenotypes[calls == g]) for g in [0, 1, 2] ] + \
               ["increases" if float(stats["BETA"]) > 0 else "decreases"]

    def carrier_values(self, phenotypes, isCarrier, hasCall):
        return [QuantitativeTrait.mean(phenotypes[isCarrier]), QuantitativeTrait.mean(phenotypes[hasCall & ~isCarrier])]

    def family_values(self, phenotypes):
        return [QuantitativeTrait.mean(phenotypes)]

    def family_carrier_values(self, phenotypes, isCarrier, hasCall):
        return self.carrier_values(phenotypes, isCarrier, hasCall)

TRAIT_TYPES = { "binary": BinaryTrait, "quantitative": QuantitativeTrait }

def get_trait(traitType, **settings):
    '''
    Returns the trait handler of a trait type, with its settings (e.g., as
    recorded by a 'saige gwas' manifest).
    '''
    if traitType not in TRAIT_TYPES:
        raise ValueError(f"Trait type '{traitType}' must be one of {list(TRAIT_TYPES)}")
    return TRAIT_TYPES[traitType](**settings)
