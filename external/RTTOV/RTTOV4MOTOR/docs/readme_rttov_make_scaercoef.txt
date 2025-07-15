!!   --------------------------------------------------------------------------
!!   Instructions for using rttov_make_scaercoef.exe
!!   --------------------------------------------------------------------------
!!
!!   This executable allows users to generate their own "scaercoef" optical
!!   property files.
!!
!!   Usage: rttov_make_scaercoef.exe --config_file     STRING
!!                                   --rtcoef_file     STRING
!!                                   --output_file     STRING
!!                                   --normalise_size           (OPTIONAL)
!!                                   --nthreads        INTEGER  (OPTIONAL)
!!                                   --chou_only                (OPTIONAL)
!!                                   --max_nmom        INTEGER  (OPTIONAL)
!!
!!   where config_file    - the path to a file which defines the aerosol
!!                          properties (see below)
!!         rtcoef_file    - the path to an RTTOV optical depth coefficient file
!!         output_file    - the output aerosol property file name
!!         normalise_size - if this argument is present the size distributions
!!                          are explicitly normalised to 1 particle/cm^3. By
!!                          default this is not done.
!!         nthreads       - the number of threads to use (optional, defaults
!!                          to one, recommended to compile with OpenMP
!!                          enabled and to make use of multiple threads)
!!         chou_only      - if specified the output file will only contain
!!                          optical properties required for IR Chou-scaling
!!                          simulations: much smaller file, but cannot be
!!                          used with the DOM solver or for solar simulations
!!                          (default false i.e. generate optical properties
!!                          for all solvers)
!!         max_nmom       - the maximum number of Legendre coefficients to
!!                          calculate for each phase function; this defines
!!                          the max number of DOM streams that can be used in
!!                          RTTOV with these optical properties, default 128
!!
!!   Optical properties are computed using Mie theory. Optical properties can
!!   be calculated for one or more particle types. Each particle type may be
!!   hydrophobic (no dependence on relative humidity) or hydrophilic. In the
!!   latter case optical properties are calculated for multiple values of
!!   relative humidity (RH): within RTTOV these are interpolated to the RH of
!!   the layer.
!!
!!   Optical properties are computed for all channels in the rtcoef_file. If
!!   any channels are solar-enabled (applies only to "v9 predictor" coef files)
!!   then phase functions will be stored for these channels unless the
!!   chou_only argument is passed. If you only need optical properties for a
!!   subset of channels then use the rttov_conv_coef.exe executable to extract
!!   these channels first. The format of the output aerosol file is the same as
!!   that of the supplied rtcoef_file (ASCII, binary or HDF5).
!!
!!   Note that the number of levels of the coefficient file and the predictor
!!   version do not affect the calculated optical properties (except in the
!!   case of solar-enabled channels in v9 predictor files). Therefore you can
!!   use the same aerosol file with v7, v8 or v9 IASI coefficients on 54L or
!!   101L for non-solar calculations. If you want to include solar radiation
!!   you would need to generate an aerosol file using a v9 predictor IASI file.
!!
!!   The inputs required for each particle type (and for each RH value for
!!   hydrophilic species) are:
!!   - complex refractive index defined on some wavelength grid covering the
!!     spectral range of the channels in the rtcoef_file
!!   - particle size distribution defined on a suitable radius grid,
!!     normalised to 1 particle per unit volume
!!   - information about the particle density (see below)
!!
!!   This information is supplied via the config_file. A sample file is
!!   provided in data/example_rttov_make_scaercoef_config.txt.
!!   This defines the following:
!!   - number of particle types
!!   - for each particle type:
!!     - particle name
!!     - the number of relative humidity values
!!     - the relative humidity values
!!     - the factor used to convert mass mixing ratio to number density
!!       (usually for RH=0%); if <=0 it is computed from the size distribution
!!       for RH=0% and the supplied mass density value
!!     - mass density in g.cm^-3; this value is only used if the unit
!!       conversion factor is <=0
!!     - for each relative humidity value:
!!       - the path to a file containing the refractive indices
!!       - the path to a file containing the size distribution
!!
!!   There is no limit on the number of particle types.
!!
!!   The particle name is a four character field used to identify the particle.
!!
!!   For hydrophobic species the number of RH values should be set to 1 and
!!   the RH value should be 0%.
!!
!!   For hydrophilic species there is no limit to the number of RH values, but
!!   more values imply larger files. RH is specified as integer percent values
!!   starting at 0%, monotonically ascending, and not exceeding 99%. As an
!!   example, the OPAC aerosols use 8 RH values: 0, 50, 70, 80, 90, 95, 98, 99.
!!
!!   The complex refractive index data should be provided for a range of
!!   wavelengths covering the central wavenumbers of the channels in the
!!   supplied rtcoef_file. A warning is printed out if channels lie beyond the
!!   range of the data: optical properties are computed using constant-value
!!   extrapolation of the refractive indices in this case. The wavelength grid
!!   does not need to be evenly spaced. Example refractive index input files
!!   are provided in the data/ directory. The imaginary part of the refractive
!!   index (k) may be positive or negative, but is always treated as if it is
!!   negative i.e. the value -ABS(k) is used in the calculations.
!!
!!   The size distribution is specified as number density values normalised to
!!   one particle per unit volume for a range of particle radius values. The
!!   radius units are microns (um) and the normalised size distribution units
!!   are um^-1. The optical properties stored in the scaercoef files are
!!   computed per [particle.cm^-3]. Inside RTTOV they are then multiplied by
!!   the particle number density in [particles.cm^-3] to compute optical
!!   depths. If the --normalise_size argument is passed to the executable then
!!   the size distributions will be explicitly normalised by the code. The
!!   radius grid does not need to be evenly spaced. Example size distribution
!!   input files are provided in the data/ directory.
!!
!!   NB In all input files, comment lines begin with ! and are optional, and
!!   file paths must be enclosed in quotes "".
!!
!!   The conversion factor is used inside RTTOV to convert input aerosol mass
!!   mixing ratio values to particle number density. It is not used if the
!!   input aerosol units are number density so if you are only going to supply
!!   aerosol concentrations in number density then the value of this conversion
!!   factor does not matter (e.g. set it to 1).
!!
!!   The units of the conversion factor are [g.m^-3]/[particle.cm^-3] and it
!!   represents the mass density per particle density. Given aerosol mass
!!   density rho, the conversion factor is calculated as:
!!
!!          rho * 4pi/3 * integral[r^3.size_dist(r).dr]
!!
!!   where the integral is over the range of radius values for which the size
!!   distribution is specified in microns.
!!
!!   The conversion factor may be calculated and specified explicitly in the
!!   config_file. For example, it may be desirable to limit the calculation to
!!   a maximum particle size: this is done for the OPAC aerosols. Alternatively
!!   if the conversion factor in the config_file is set to a value less than or
!!   equal to zero then a value for the density (rho) of the aerosol material
!!   must be supplied in units of [g.cm^-3]. The conversion factor is then
!!   computed from the size distribution for RH=0% and this density value.
!!
!!   NB The value of confac stored in the scaercoef file is the reciprocal of
!!   the value given by the expression above, or, if greater than zero, the
!!   reciprocal of the value explicitly specified in the config_file.
!!
!!   You can test this executable using the example input data from the
!!   top-level RTTOV directory as follows:
!!
!!   $ bin/rttov_make_scaercoef.exe \
!!   > --config_file data/example_rttov_make_scaercoef_config.txt \
!!   > --rtcoef_file rtcoef_rttov13/rttov7pred54L/rtcoef_msg_3_seviri.dat \
!!   > --output_file scaercoef_msg_3_seviri_example.dat \
!!   > --normalise_size \
!!   > --nthreads 8 \
!!   > --chou_only
!!
!!   You can then compare your output file with the reference one:
!!     data/scaercoef_msg_3_seviri_example.dat
