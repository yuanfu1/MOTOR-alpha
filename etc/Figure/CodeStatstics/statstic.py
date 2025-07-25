import plotly.graph_objects as go
color_discrete_sequence = ['', '#1abc9c', '#2ecc71', '#3498db', '#9b59b6',
                           '#34495e', '#f1c40f', '#e67e22', '#e74c3c', '#ecf0f1', '#95a5a6']

fig = go.Figure(go.Sunburst(
    ids=[
        "MOTOR", "MOTOR-DA", "MOTOR-QC", "MOTOR-Rep", "Utilities", "ThirdParty",
        "BMat", "IO", "MGOpts", "Minimization", "Obs", "RMat", "States", "Transform",
        "ObsBase", "ObsConvention", "ObsRadar", "ObsSatellite", "ObsUtilities",
        "AuxType", "Geometry", "MPDD", "MPObs", "StructureGrid",
        "FLog", "GrapesIO", "InterpHP", "NMLRead", "Utility",
        "LBFGSB", "BLAS", "LinPACK", "LAPACK", "NetCDF", "OpenMPI",
        "BField", "BMatrix", "EnLoc",
        "IOModels", "IOMOTOR",
        "JFunc", "MinSolver",
        "C2O", "ObsField", "ObsSet",
        "RMatrix", "RField",
        "Field", "State",
        "Constraint", "CVOperator", "M2ODirect", "ObsOperator", "TransBase",
        "AdvanceTime", "Relaxations", "WRFIO", "interp1D",
        "ObsSound", "ObsSurface", "ObsVwpw",
        "ObsRadarRef", "ObsRadarVel",
        "MultiGrid", "SingleGrid",
        "AuxTypeObs", "AuxTypeSg",
        "Application",
        "3DVarGrapes", "3DVarWRF", "3DVarGrapesMG", "3DVarWRFMG",
        "UniTest",
        "Convertion","SFCRandoms"
    ],
    labels=[
        "MOTOR", "MOTOR-DA", "MOTOR-QC", "MOTOR-Rep", "Utilities", "ThirdParty",
        "BMat", "IO", "MGOpts", "Minimization", "Obs", "RMat", "States", "Transform",
        "ObsBase", "ObsConvention", "ObsRadar", "ObsSatellite", "ObsUtilities",
        "AuxType", "Geometry", "MPDD", "MPObs", "StructureGrid",
        "FLog", "GrapesIO", "InterpHP", "NMLRead", "Utility",
        "LBFGSB", "BLAS", "LinPACK", "LAPACK", "NetCDF", "OpenMPI",
        "BField", "BMatrix", "EnLoc",
        "IOModels", "IOMOTOR",
        "JFunc", "MinSolver",
        "C2O", "ObsField", "ObsSet",
        "RMatrix", "RField",
        "Field", "State",
        "Constraint", "CVOperator", "M2ODirect", "ObsOperator", "TransBase",
        "AdvanceTime", "Relaxations", "WRFIO", "interp1D",
        "ObsSound", "ObsSurface", "ObsVwpw",
        "ObsRadarRef", "ObsRadarVel",
        "MultiGrid", "SingleGrid",
        "AuxTypeObs", "AuxTypeSg",
                "Application",
        "3DVarGrapes", "3DVarWRF", "3DVarGrapesMG", "3DVarWRFMG",
                "UniTest",
        "Convertion", "SFCRandoms"

    ],
    parents=[
        "", "MOTOR", "MOTOR", "MOTOR", "MOTOR", "MOTOR",
        "MOTOR-DA", "MOTOR-DA", "MOTOR-DA", "MOTOR-DA", "MOTOR-DA", "MOTOR-DA", "MOTOR-DA", "MOTOR-DA",
        "MOTOR-QC", "MOTOR-QC", "MOTOR-QC", "MOTOR-QC", "MOTOR-QC",
        "MOTOR-Rep", "MOTOR-Rep", "MOTOR-Rep", "MOTOR-Rep", "MOTOR-Rep",
        "Utilities", "Utilities", "Utilities", "Utilities", "Utilities",
        "ThirdParty", "ThirdParty", "ThirdParty", "ThirdParty", "ThirdParty", "ThirdParty",
        "BMat", "BMat", "BMat",
        "IO", "IO",
        "Minimization", "Minimization",
        "Obs", "Obs", "Obs",
        "RMat", "RMat",
        "States", "States",
        "Transform", "Transform", "Transform", "Transform", "Transform",
        "Utility", "Utility", "Utility", "Utility",
        "ObsConvention", "ObsConvention", "ObsConvention",
        "ObsRadar", "ObsRadar",
        "Geometry", "Geometry",
        "AuxType", "AuxType",
                "MOTOR",
                        "Application", "Application", "Application", "Application",
                        "Application",
        "UniTest", "UniTest"
    ], marker=dict(colors=color_discrete_sequence)
))
fig.update_layout(margin=dict(t=0, l=0, r=0, b=0))
fig.show()
