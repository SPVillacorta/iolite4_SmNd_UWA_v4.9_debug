#/ Type: DRS
#/ Name: Sm-Nd Isotopes_UWA_v4.8
#/ Authors: S Villacorta, C Fisher
#/ Description: Sm-Nd DRS script for Iolite 4. adds visible rho (AssociatedResult + time series), spline validation, and error propagation.
#/ References: Fisher et al. 2020; Barrote et al. 2021
#/ Version: 4.8
#/ Contact: villacortasp@gmail.com


from iolite import QtGui
from iolite.Qt import Qt, QColor
from iolite.ui import IolitePlotPyInterface as Plot
from iolite.types import Result
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import numpy as np


def rho_corr_147Sm_143Nd(sel):
    """
    Calculates the Pearson correlation coefficient (rho) between log10 of
    standard-corrected 147Sm/144Nd and 143Nd/144Nd ratios for a given selection.
    """
    result = Result()
    try:
        sm = data.timeSeries("StdCorr_147Sm_144").dataForSelection(sel)
        nd = data.timeSeries("StdCorr_143Nd_144").dataForSelection(sel)
        # Ensure only valid, non-zero data points are used for correlation
        valid = ~np.isnan(sm) & ~np.isnan(nd) & (sm > 0) & (nd > 0)
        if np.any(valid):
            rho_val = np.corrcoef(np.log10(sm[valid]), np.log10(nd[valid]))[0, 1]
            result.setValue(rho_val)
        else:
            result.setValue(np.nan)
    except Exception as e:
        drs.message(f"[rho_corr_147Sm_143Nd] Error in selection {sel.name()}: {e}")
        result.setValue(np.nan)
    return result

def runDRS():
    drs.message("Running Sm-Nd Isotopes_UWA_v48...")
    drs.progress(0)

    settings = drs.settings()
    indexChannel = data.timeSeries(settings["IndexChannel"])
    time = indexChannel.time()
    maskChannel = data.timeSeries(settings["MaskChannel"])
    rmName = settings["ReferenceMaterial"]
    rmNames = data.selectionGroupNames(data.ReferenceMaterial)

    # check for Reference Material
    if rmName not in rmNames:
        if "G_LREE" in rmNames:
            rmName = "G_LREE"
        elif rmNames:
            rmName = rmNames[0]
        else:
            drs.message("Error: No reference material available. Please define one in your data.")
            drs.finished()
            return # Exit if no RM is found
        drs.message(f"Selected reference material '{settings['ReferenceMaterial']}' not found. Using '{rmName}' instead.")
        drs.setSetting("ReferenceMaterial", rmName) # Update the setting to the chosen RM

    cutoff = settings["MaskCutoff"]
    trim = settings["MaskTrim"]
    NdTrue = settings["NdTrue"]
    Sm147_149 = settings["Sm147_149"]
    Sm144_149 = settings["Sm144_149"]
    Sm148_149 = settings["Sm148_149"]
    propErrors = settings["PropagateError"]

    drs.setIndexChannel(indexChannel)
    mask = drs.createMaskFromCutoff(maskChannel, cutoff, trim)

    # Baseline Subtraction
    for ch in data.timeSeriesList(data.Input):
        drs.baselineSubtract(data.selectionGroup("Baseline"), [ch], mask, 25, 75)

    # Retrieve CPS data
    Nd143 = data.timeSeries("Nd143_CPS").data()
    Nd144 = data.timeSeries("Nd144_CPS").data() # Note: This will be corrected for Sm contribution later
    Nd145 = data.timeSeries("Nd145_CPS").data()
    Nd146 = data.timeSeries("Nd146_CPS").data()
    Nd148 = data.timeSeries("Nd148_CPS").data() # Note: This will be corrected for Sm contribution later
    Sm147 = data.timeSeries("Sm147_CPS").data()
    Sm149 = data.timeSeries("Sm149_CPS").data()

    # Sm Mass Fractionation Correction
    SmFract = (np.log(Sm147_149 / np.where(Sm149 != 0, Sm147 / Sm149, np.nan))) / np.log(146.9149 / 148.9172) * mask
    Sm144 = Sm149 * Sm144_149 * ((148.9172 / 143.912) ** SmFract) * mask
    Sm148 = Sm149 * Sm148_149 * ((148.9172 / 147.9148) ** SmFract) * mask

    # Nd interference correction
    Nd144c = (Nd144 - Sm144) * mask
    Nd148c = (Nd148 - Sm148) * mask
    # Ensure Nd144c is not zero for subsequent divisions
    Nd144c_safe = np.where(Nd144c != 0, Nd144c, np.nan)

    # Nd Mass Fractionation Correction
    # Robust division by Nd144c_safe for NdFract calculation
    NdFract = (np.log(NdTrue / np.where(Nd144c_safe != 0, Nd146 / Nd144c_safe, np.nan))) / np.log(145.9131 / 143.9098) * mask

    # Calculate corrected ratios
    Nd143_144_Corr = (Nd143 / Nd144c_safe) * (142.9098 / 143.9098) ** NdFract * mask
    Nd145_144_Corr = (Nd145 / Nd144c_safe) * (144.9126 / 143.9098) ** NdFract * mask
    Nd148_144_Corr = (Nd148c / Nd144c_safe) * (147.916889 / 143.9098) ** NdFract * mask
    Sm147_Nd144_Corr = (Sm147 / Nd144c_safe) * mask # Raw Sm147/Nd144 ratio after Nd correction

    # Total Nd calculation
    TotalNd = Nd146 / 0.172 * mask

    # Create uncorrected time series for spline generation
    data.createTimeSeries("Nd143_144_Corr", data.Output, time, Nd143_144_Corr)
    data.createTimeSeries("Sm147_Nd144_Corr", data.Output, time, Sm147_Nd144_Corr)

    # Spline-based standard correction for 143Nd/144Nd
    StdCorr_143Nd_144 = np.full_like(Nd143_144_Corr, np.nan) # Initialize with NaN
    try:
        data.createSpline(rmName, "Nd143_144_Corr")
        spline143 = data.spline(rmName, "Nd143_144_Corr").data()
        std143 = data.referenceMaterialData(rmName)["143Nd/144Nd"].value()
        StdCorr_143Nd_144 = Nd143_144_Corr * std143 / np.where(spline143 != 0, spline143, np.nan)
    except Exception as e:
        drs.message(f"Warning: Standard correction for 143Nd/144Nd failed: {e}")

    # Spline-based standard correction for 147Sm/144Nd (calculation logic as 143Nd/144Nd)
    StdCorr_147Sm_144 = np.full_like(Sm147_Nd144_Corr, np.nan) # Initialize with NaN
    try:
        data.createSpline(rmName, "Sm147_Nd144_Corr")
        spline147 = data.spline(rmName, "Sm147_Nd144_Corr").data()
        std147 = data.referenceMaterialData(rmName)["147Sm/144Nd"].value()
        StdCorr_147Sm_144 = Sm147_Nd144_Corr * std147 / np.where(spline147 != 0, spline147, np.nan)
    except Exception as e:
        drs.message(f"Warning: Standard correction for 147Sm/144Nd failed: {e}")

    # Create Output Time Series (in specified order for publication)
    data.createTimeSeries("StdCorr_147Sm_144", data.Output, time, StdCorr_147Sm_144)
    data.createTimeSeries("StdCorr_143Nd_144", data.Output, time, StdCorr_143Nd_144)
    data.createTimeSeries("Nd143_144_Corr", data.Output, time, Nd143_144_Corr)
    data.createTimeSeries("Nd145_144_Corr", data.Output, time, Nd145_144_Corr)
    data.createTimeSeries("Nd148_144_Corr", data.Output, time, Nd148_144_Corr)
    data.createTimeSeries("NdFract", data.Output, time, NdFract)
    data.createTimeSeries("SmFract", data.Output, time, SmFract)
    data.createTimeSeries("TotalNd", data.Output, time, TotalNd)

    # Error Propagation
    if propErrors:
        selGroups = [g for g in data.selectionGroupList() if g.type != data.Baseline]
        data.propagateErrors(selGroups, [data.timeSeries("StdCorr_143Nd_144")], data.timeSeries("Nd143_144_Corr"), rmName)
        # Propagate errors for StdCorr_147Sm_144 as well
        data.propagateErrors(selGroups, [data.timeSeries("StdCorr_147Sm_144")], data.timeSeries("Sm147_Nd144_Corr"), rmName)
        drs.message("Propagated uncertainties calculated.")

    # Register Associated Result for Rho
    # rho from StdCorr_147Sm_Nd144 and StdCorr_143Nd_144
    data.registerAssociatedResult("log10(StdCorr_147Sm/144Nd) vs log10(StdCorr_143Nd/144Nd) rho", rho_corr_147Sm_143Nd)

    # To ensure 'rho' is always visible in the output panel
    # For publication, the primary rho value will be in the associated results,
    # but a flat line helps in plotting/visibility if desired.
    rho_ts_data = np.full_like(time, np.nan)
    data.createTimeSeries("rho", data.Output, time, rho_ts_data)
    # Update this time series with the first valid rho found
    first_rho_val = np.nan
    for sel in data.selectionGroupList():
        if sel.type != data.Baseline:
            try:
                v = rho_corr_147Sm_143Nd(sel).value()
                if not np.isnan(v):
                    first_rho_val = v
                    break
            except:
                continue
    if not np.isnan(first_rho_val):
        data.timeSeries("rho").setData(np.full_like(time, first_rho_val))

    drs.progress(100)
    drs.message("Finished!")
    drs.finished()

def settingsWidget():
    widget = QtGui.QWidget()
    formLayout = QtGui.QFormLayout()
    formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
    formLayout.setFormAlignment(Qt.AlignHCenter | Qt.AlignTop)
    widget.setLayout(formLayout)

    timeSeriesNames = data.timeSeriesNames(data.Input)
    rmNames = data.selectionGroupNames(data.ReferenceMaterial)

    # Default settings (using v1.35 defaults where appropriate)
    drs.setSetting("IndexChannel", "Nd146")
    drs.setSetting("ReferenceMaterial", "G_LREE" if "G_LREE" in rmNames else (rmNames[0] if rmNames else ""))
    drs.setSetting("MaskChannel", "Nd146")
    drs.setSetting("MaskCutoff", 0.1)
    drs.setSetting("MaskTrim", 0.0) # Retain MaskTrim setting
    drs.setSetting("NdTrue", 0.7219)
    drs.setSetting("Sm147_149", 1.08680)
    drs.setSetting("Sm144_149", 0.22332)
    drs.setSetting("Sm148_149", 0.81419)
    drs.setSetting("PropagateError", True)

    settings = drs.settings()

    def addLine(label, val, key):
        le = QtGui.QLineEdit(widget)
        le.setText(str(val))
        # Use a validator to ensure float input
        validator = QtGui.QDoubleValidator()
        le.setValidator(validator)
        le.textChanged.connect(lambda t: drs.setSetting(key, float(t) if t else 0.0)) # Handle empty string
        formLayout.addRow(label, le)

    def addCombo(label, items, val, key):
        cb = QtGui.QComboBox(widget)
        cb.addItems(items)
        cb.setCurrentText(val)
        cb.currentTextChanged.connect(lambda v: drs.setSetting(key, v))
        formLayout.addRow(label, cb)

    # Populate UI fields
    addCombo("Index channel", timeSeriesNames, settings["IndexChannel"], "IndexChannel")
    addCombo("Reference material", rmNames, settings["ReferenceMaterial"], "ReferenceMaterial")
    addCombo("Mask channel", timeSeriesNames, settings["MaskChannel"], "MaskChannel")

    addLine("Mask cutoff", settings["MaskCutoff"], "MaskCutoff")
    addLine("Mask trim", settings["MaskTrim"], "MaskTrim") # Re-added for UI visibility
    addLine("NdTrue (146Nd/144Nd)", settings["NdTrue"], "NdTrue")
    addLine("147Sm/149Sm Ratio", settings["Sm147_149"], "Sm147_149")
    addLine("144Sm/149Sm Ratio", settings["Sm144_149"], "Sm144_149")
    addLine("148Sm/149Sm Ratio", settings["Sm148_149"], "Sm148_149")

    propagateBox = QtGui.QCheckBox(widget)
    propagateBox.setChecked(settings["PropagateError"])
    propagateBox.toggled.connect(lambda t: drs.setSetting("PropagateError", bool(t)))
    formLayout.addRow("Propagate Error", propagateBox)

    drs.setSettingsWidget(widget)