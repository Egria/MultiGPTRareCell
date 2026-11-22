#!/bin/bash
set -e
project=$1
r_module=$2
python=$3
shift 3

if ! type module >/dev/null 2>&1; then
    for init in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
        if [[ -r "$init" ]]; then
            source "$init" >&2
            if type module >/dev/null 2>&1; then break; fi
        fi
    done
fi
if ! type module >/dev/null 2>&1; then
    echo "RareCell: cannot initialize environment modules to load $r_module" >&2
    exit 2
fi
module load "$r_module" >&2
if ! rscript=$(type -P Rscript); then
    # The DCC RHEL8 module is also used from newer login-node images.
    if [[ "$r_module" == "R/4.3.1-rhel8" && -x /opt/apps/rhel8/R-4.3.1/bin/Rscript ]]; then
        export PATH="/opt/apps/rhel8/R-4.3.1/bin:$PATH"
        rscript=/opt/apps/rhel8/R-4.3.1/bin/Rscript
    else
        echo "RareCell: module $r_module did not expose an executable Rscript" >&2
        exit 2
    fi
fi
[[ "$rscript" = /* ]] || rscript="$PWD/$rscript"
export RARECELL_MODULE_RSCRIPT="$rscript"
export R_LIBS_USER="${RARECELL_R_LIBS_USER:-$(dirname "$project")/r-libs/MultiGPTRareCell/4.3.1-rhel8}"
exec "$python" -m rarecell.cli "$@"
