#!/bin/bash
if [[ "$(uname -s)" == "Linux" || "$(uname -s)" == "Darwin" ]]; then
    export PATH="$HOME/bin/sratoolkit/bin:$PATH"
    vdb-config --set "/repository/user/main/public/root=$HOME/.ncbi/public"
else
    export PATH="$(cygpath -w "$HOME/bin/sratoolkit/bin"):$PATH"
    vdb-config --set "/repository/user/main/public/root=$(cygpath -w "$HOME/.ncbi/public")"
fi

