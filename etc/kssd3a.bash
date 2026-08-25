# bash completion for kssd3a

_kssd3a()
{
    local cur prev cmd opts

    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    if [[ ${COMP_CWORD} -eq 1 ]]; then
        opts="sketch ani matrix set examples doctor help shuffle dist reverse --help --usage --version --license"
        COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
        return 0
    fi

    cmd="${COMP_WORDS[1]}"

    case "${cmd}" in
        sketch)
            case "${prev}" in
                -o|--outdir|-l|--list|-i|--index|--remove|--keep)
                    compopt -o filenames 2>/dev/null
                    COMPREPLY=( $(compgen -f -- "${cur}") )
                    return 0
                    ;;
                --metric)
                    opts="ctx-moe ctx-naive mash aaf ctx-naive&aaf ctx-moe|mash"
                    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
                    return 0
                    ;;
            esac
            opts="-l --list -o --outdir -P --pipecmd -p --threads -T --use_coden_ctxobj -C --ctxlen -O --outerobjlen -I --innerobjlen -f --DimRdcFold -n --LstKmerOcrs --npercentile --ncap --readsQC --sketchQC -A --abundance --anno --position --conflict -a --asone --splitmfa --psmp --psketch --pindex --ppos -i --index --append --remove --keep --drop-position --dedup --dedup-max-afcut --dedup-ctxcut --dedup-index --dedup-index-max-ctx-freq --dedup-index-min-votes --dedup-index-sample-step --metric --merge --help --usage"
            ;;
        ani)
            case "${prev}" in
                -r|--ref|-q|--query|--qraw|--reflist|--qrylist|--pipecmd|-g|--glist|-o|--outfile)
                    compopt -o filenames 2>/dev/null
                    COMPREPLY=( $(compgen -f -- "${cur}") )
                    return 0
                    ;;
            esac
            opts="-r --ref -q --query --qraw --reflist --qrylist -v --naive -m --outfmt -s --slmetrics --unified-metric -t --ctxcut -f --afcut -n --anicut -c --control -N --top -g --glist -o --outfile -p --threads -d --diagonal -e --exception --pair --use_coden_ctxobj --DimRdcFold --ctxlen --outerobjlen --innerobjlen --LstKmerOcrs --npercentile --ncap --readsQC --abundance --anno --conflict --asone --splitmfa --pipecmd --help --usage"
            ;;
        matrix)
            case "${prev}" in
                -r|--ref|-q|--query|-g|--glist|-o|--outfile|--edge-out|--keep-out|--remove-out|--matrix-idmap|--keep-matrix-out|--keep-matrix-idmap)
                    compopt -o filenames 2>/dev/null
                    COMPREPLY=( $(compgen -f -- "${cur}") )
                    return 0
                    ;;
                --matrix-format|--keep-matrix-format)
                    opts="tsv phylip"
                    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
                    return 0
                    ;;
                --progress)
                    opts="auto on off"
                    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
                    return 0
                    ;;
                -m|--metric)
                    opts="ctx-moe ctx-naive mash aaf 0 1 ctx-naive&aaf ctx-moe|mash"
                    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
                    return 0
                    ;;
                --format)
                    opts="full triangle edges clusters dedup-plan"
                    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
                    return 0
                    ;;
            esac
            opts="-r --ref -q --query -m --metric --format --cut --max-afcut --ctxcut --index-max-ctx-freq --index-min-votes --index-sample-step --progress -g --glist -o --outfile --matrix-format --matrix-idmap --edge-out --keep-out --remove-out --keep-matrix-out --keep-matrix-format --keep-matrix-idmap -p --threads -d --diagonal --diagonal-value -e --exception --help --usage"
            ;;
        set)
            case "${prev}" in
                -o|--outdir|-g|--grouping|-s|--subtract|-i|--intsect)
                    compopt -o filenames 2>/dev/null
                    COMPREPLY=( $(compgen -f -- "${cur}") )
                    return 0
                    ;;
            esac
            opts="-u --union -s --subtract -i --intsect -q --uniq_union --markerdb -p --threads -P --print --psketch --pindex --ppos -g --grouping -o --outdir --help --usage"
            ;;
        dist)
            case "${prev}" in
                -l|--list|-r|--reference_dir|-o|--outdir|-f|--skf)
                    compopt -o filenames 2>/dev/null
                    COMPREPLY=( $(compgen -f -- "${cur}") )
                    return 0
                    ;;
            esac
            opts="-k --halfKmerlength -p --threadN -l --list -L --DimRdcLevel -m --maxMemory -n --LstKmerOcrs -Q --quality -r --reference_dir -o --outdir -N --neighborN_max -D --mutDist_max -M --metric -C --containment -O --outfields --correction -A --abundance -u --dedup --keepcofile -P --pipecmd --keepskf -f --skf --byread --help --usage"
            ;;
        *)
            opts="--help --usage --license"
            ;;
    esac

    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
    return 0
}

complete -F _kssd3a kssd3a
