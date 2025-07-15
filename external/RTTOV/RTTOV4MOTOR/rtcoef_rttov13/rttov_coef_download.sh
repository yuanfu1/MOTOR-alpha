#!/bin/bash

# =================================================================

# Script to download RTTOV coefficient files from the website:
# run the script from within your rtcoef_rttov13/ directory.

# The function "download_extract" downloads all files with a given
# file extension from the specified directory. It then unzips and
# untars the files as necessary. Intermediate (e.g. .tar) files
# are deleted so only coefficient files are left.

# It uses the local directory structure of rtcoef_rttov13/ to
# determine the rttov*pred*L/ directories on the server.

# Filetypes are distinguished as follows:
# - HDF5 files end in .H5
# - tarballs of smaller coef files end in .tar.bz2
# - other large inidividually zipped files end in .dat.bz2

# =================================================================


echo '========================================================================='
echo ' RTTOV coefficient download script'
echo '========================================================================='
echo
echo '* Individual coefficient files may be downloaded from the RTTOV website:'
echo '  http://nwpsaf.eu/site/software/rttov/downloads/coefficients/coefficient-download/'
echo '  See the website for information about the different types of coefficient'
echo '  files.'
echo
echo '* If you only require specific coefficient files it may be quicker to'
echo '  download them from the website directly.'
echo
echo '* This script downloads all coefficient files or all files of a particular'
echo '  type. You must call this from within the rtcoef_rttov13 directory:'
echo '  the coefficient files are unzipped and extracted in the standard'
echo '  locations within this directory.'
echo
echo '* Be aware: EXISTING COEFFICIENT FILES MAY BE OVERWRITTEN WITHOUT WARNING.'
echo
echo '* You only need to download the coefficient files for the instruments you'
echo '  are simulating and the kinds of simulations you are carrying out.'
echo
echo '* The v13 predictor VIS/IR/MW optical depth coefficient files are included'
echo '  with the RTTOV distribution. However the tarballs on the website may'
echo '  contain updated coefficient files. This script downloads all available'
echo '  v13 predictor and v7/v8/v9 predictor rtcoef files.'
echo
echo ' * If a download fails part-way through you can request that the script'
echo '   attempts to continue the download rather than starting again.'
echo
echo ' * In order to enable continuation of failed downloads the script does'
echo '   not delete the temporary download directory. You can remove this once'
echo '   your downloads are complete.'
echo
echo 'File types are:'
echo '  - Optical depth (rtcoef) coefficient files for MW, IR and VIS/IR'
echo '  - VIS/IR aerosol (scaer) and cloud (sccld) coefficient files'
echo '  - MFASIS LUT files'
echo '  - RTTOV-SCATT MW hydrotable files'
echo '  - Optical depth (rtcoef) files for hi-res IR sounders'
echo '  - Hi-res IR sounder aerosol (scaer) and cloud (sccld) coefficient files'
echo '  - PC-RTTOV coefficient files'
echo '  - HTFRTC coefficient files in netCDF format'
echo '  - HTFRTC coefficient files in ASCII format'
echo
echo '========================================================================='
echo


verbose_wget=0 # Set this to a non-zero value to see full wget output

tmpdir='tmp_downloads/'
protocol='https://'
server='nwp-saf.eumetsat.int/'
path='downloads/rtcoef_rttov13/'

download_extract () {
    # Downloads all files from a directory on server with given file extension
    # Zipped files are unzipped, tarballs are extracted and .tar files deleted
    # Tarball extraction overwrites existing files
    # Other files will be renamed rather than being overwritten

    # First argument is file extension e.g. ".tar.bz2" or ".dat.bz2" or ".H5" or ".nc"
    # Second argument is directory e.g. "cldaer_visir" or "rttov13pred54L"

    d=$2
    echo
    echo "=============================================================================="
    echo "Downloading files with extension $1 from ${protocol}${server}${path}${d}"
    echo
    if [ $verbose_wget -eq 0 ]; then
        wget "${continue}" -q -r -np -l1 -A"$1" ${protocol}${server}${path}"${d}" -P"${tmpdir}"
    else
        echo "--- wget output start ------------------------------------------------------"
        echo
        wget "${continue}" -r -np -l1 -A"$1" ${protocol}${server}${path}"${d}" -P"${tmpdir}"
        echo
        echo "--- wget output end --------------------------------------------------------"
        echo
    fi
    local result=$?
    if [ $result -ne 0 ]; then
        echo "wget failure with error code $result"
        return
    fi
    if [[ -d ${tmpdir}${server}${path}"${d}" ]]; then
        if [[ $(ls ${tmpdir}${server}${path}"${d}") ]]; then
            for f in ${tmpdir}${server}${path}"${d}"/*; do
                echo "File found: $(basename "$f")"
                if [[ -f ${d}/$(basename "$f") ]]; then
                    mv "${d}/$(basename "$f")" "${d}/$(basename "$f")_old_$(date +%Y%m%d_%H%M%S)"
                    echo "Renamed existing file ${d}/$(basename "$f")"
                fi
            done
            mv ${tmpdir}${server}${path}"${d}"/* "$d"
        else
            echo "No files found with names: ${protocol}${server}${path}${d}/*$1"
            echo "=============================================================================="
            echo
            return
        fi
    else
        echo "Directory not found on server: ${protocol}${server}${path}${d}"
        echo "=============================================================================="
        echo
        return
    fi

    echo
    echo "Extracting downloaded files..."
    echo

    cd "$d"
    for f in *; do
        if [[ $f = *".bz2" ]]; then
            echo "Unzipping $f..."
            bunzip2 "$f" 2> /dev/null
            local result=$?
            if [ $result -eq 1 ]; then
                echo "Unzipped file is already present, not unzipping $f"
            elif [ $result -ne 0 ]; then
                echo "Error unzipping, not unzipping $f"
            fi
        fi
    done

    for f in *; do
        if [[ $f = *".tar" ]]; then
            echo "Extracting tarball $f..."
            tar xf "$f" 2> /dev/null
            local result=$?
            if [ $result -ne 0 ]; then
                echo "Tar extraction failed, not extracting $f"
            fi
            rm -f "$f"
        fi
    done
    echo "=============================================================================="
    echo
    cd ../
}


# Check we are in rtcoef_rttov13/ directory; if not then try to find it

if [[ $(pwd) != *"rtcoef_rttov13" ]]; then
    if [[ -d rtcoef_rttov13 ]]; then
        echo "Not in RTTOV rtcoef_rttov13/ directory, changing to rtcoef_rttov13/"
        cd rtcoef_rttov13   
    elif [[ -d ../rtcoef_rttov13 ]]; then
        echo "Not in RTTOV rtcoef_rttov13/ directory, changing to rtcoef_rttov13/"
        cd ../rtcoef_rttov13
    elif [[ -d ../../rtcoef_rttov13 ]]; then
        echo "Not in RTTOV rtcoef_rttov13/ directory, changing to rtcoef_rttov13/"
        cd ../../rtcoef_rttov13
    else
        echo "This script must be run from the RTTOV rtcoef_rttov13/ directory."
        exit 1
    fi
    echo
fi

# echo "For the larger coefficient files do you want HDF5 (y) or ASCII (n) format? "
# while true; do
#     read -p "> " yn
#     case $yn in
#         [Yy] ) hdf=1; ascii=0; break;;
#         [Nn] ) hdf=0; ascii=1; break;;
#         * ) echo "Please answer y or n.";;
#     esac
# done
hdf=1
ascii=0


echo "Download all files (y) or specify files to download (n)? "
while true; do
    read -p "> " yn
    case $yn in
        [Yy] ) yn=1; break;;
        [Nn] ) yn=0; break;;
        * ) echo "Please answer y or n.";;
    esac
done
if [ $yn -eq 1 ]; then
    get_visirmw=1
    get_ircldaer=1
    get_mfasis=1
    get_hydrotables=1
    get_hiresascii=$ascii
    get_hireshdf5=$hdf
    get_hirescldaerascii=$ascii
    get_hirescldaerhdf5=$hdf
    get_pcascii=$ascii
    get_pchdf5=$hdf
    get_htfrtc_netcdf=1
    get_htfrtc_ascii=1
else
    echo "Download VIS/IR/MW rtcoef files? (y/n) "
    while true; do
        read -p "> " yn
        case $yn in
            [Yy] ) get_visirmw=1; break;;
            [Nn] ) get_visirmw=0; break;;
            * ) echo "Please answer y or n.";;
        esac
    done

    echo "Download VIS/IR cld/aer coef files? (y/n) "
    while true; do
        read -p "> " yn
        case $yn in
            [Yy] ) get_ircldaer=1; break;;
            [Nn] ) get_ircldaer=0; break;;
            * ) echo "Please answer y or n.";;
        esac
    done

    echo "Download MFASIS LUT files? (y/n) "
    while true; do
        read -p "> " yn
        case $yn in
            [Yy] ) get_mfasis=1; break;;
            [Nn] ) get_mfasis=0; break;;
            * ) echo "Please answer y or n.";;
        esac
    done

    echo "Download MW hydrotable files? (y/n) "
    while true; do
        read -p "> " yn
        case $yn in
            [Yy] ) get_hydrotables=1; break;;
            [Nn] ) get_hydrotables=0; break;;
            * ) echo "Please answer y or n.";;
        esac
    done

    echo "Download hi-res IR rtcoef files? (y/n) "
    while true; do
        read -p "> " yn
        case $yn in
            [Yy] ) get_hireshdf5=$hdf; get_hiresascii=$ascii; break;;
            [Nn] ) get_hireshdf5=0;    get_hiresascii=0; break;;
            * ) echo "Please answer y or n.";;
        esac
    done

    echo "Download hi-res IR cld/aer coef files? (y/n) "
    while true; do
        read -p "> " yn
        case $yn in
            [Yy] ) get_hirescldaerhdf5=$hdf; get_hirescldaerascii=$ascii; break;;
            [Nn] ) get_hirescldaerhdf5=0;    get_hirescldaerascii=0; break;;
            * ) echo "Please answer y or n.";;
        esac
    done

    echo "Download PC-RTTOV coef files? (y/n) "
    while true; do
        read -p "> " yn
        case $yn in
            [Yy] ) get_pchdf5=$hdf; get_pcascii=$ascii; break;;
            [Nn] ) get_pchdf5=0;    get_pcascii=0; break;;
            * ) echo "Please answer y or n.";;
        esac
    done

    echo "Download HTFRTC netCDF coef files? (y/n) "
    while true; do
        read -p "> " yn
        case $yn in
            [Yy] ) get_htfrtc_netcdf=1; break;;
            [Nn] ) get_htfrtc_netcdf=0; break;;
            * ) echo "Please answer y or n.";;
        esac
    done

    echo "Download HTFRTC ASCII coef files? (y/n) "
    while true; do
        read -p "> " yn
        case $yn in
            [Yy] ) get_htfrtc_ascii=1; break;;
            [Nn] ) get_htfrtc_ascii=0; break;;
            * ) echo "Please answer y or n.";;
        esac
    done
fi

echo "Continue previously interrupted downloads (y) or download entire files (n)? "
while true; do
    read -p "> " yn
    case $yn in
        [Yy] ) yn=1; break;;
        [Nn] ) yn=0; break;;
        * ) echo "Please answer y or n.";;
    esac
done
if [ $yn -eq 1 ]; then
  continue="--continue"
else
  continue=
  rm -fr ${tmpdir}
fi
echo

# Download VIS/IR/MW coef files:
if [ $get_visirmw -eq 1 ]; then
    for d in *; do
        if [[ -d $d && $d = "rttov"*"pred"*"L" ]]; then
            download_extract ".tar.bz2" "$d"
        fi
    done
fi

# Download non-hires VIS/IR cld/aer files:
if [ $get_ircldaer -eq 1 ]; then
    download_extract ".tar.bz2" "cldaer_ir"
    download_extract ".tar.bz2" "cldaer_visir"
fi

# Download MFASIS LUT files:
if [ $get_mfasis -eq 1 ]; then
    download_extract ".H5" "mfasis_lut"
fi

# Download hydrotable files:
if [ $get_hydrotables -eq 1 ]; then
    download_extract ".dat.bz2" "hydrotable"
fi

# Download ASCII hi-res IR sounder files:
if [ $get_hiresascii -eq 1 ]; then
    for d in *; do
        if [[ -d $d && $d = "rttov"*"pred"*"L" ]]; then
            download_extract ".dat.bz2" "$d"
        fi
    done
fi

# Download HDF5 hi-res IR sounder files:
if [ $get_hireshdf5 -eq 1 ]; then
    for d in *; do
        if [[ -d $d && $d = "rttov"*"pred"*"L" ]]; then
            download_extract ".H5" "$d"
        fi
    done
fi

# Download ASCII hires IR cld/aer files:
if [ $get_hirescldaerascii -eq 1 ]; then
    download_extract ".dat.bz2" "cldaer_ir"
    download_extract ".dat.bz2" "cldaer_visir"
fi

# Download HDF5 hires IR cld/aer files:
if [ $get_hirescldaerhdf5 -eq 1 ]; then
    download_extract ".H5" "cldaer_ir"
    download_extract ".H5" "cldaer_visir"
fi

# Download ASCII PC-RTTOV coef files:
if [ $get_pcascii -eq 1 ]; then
    download_extract ".dat.bz2" "pc"
fi

# Download HDF5 PC-RTTOV coef files:
if [ $get_pchdf5 -eq 1 ]; then
    download_extract ".H5" "pc"
fi

# Download HTFRTC netCDF coef files:
if [ $get_htfrtc_netcdf -eq 1 ]; then
    download_extract ".nc" "htfrtc"
fi

# Download HTFRTC ASCII coef files:
if [ $get_htfrtc_ascii -eq 1 ]; then
    download_extract ".dat.bz2" "htfrtc"
fi

echo
echo 'Please visit http://nwpsaf.eu/site/software/rttov/downloads/coefficients/coefficient-download/ for information about the coefficient files'
echo
