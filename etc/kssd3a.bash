# Bash completion; option names come from the installed executable's help.
_kssd3a()
{
    local cur prev cmd opts
    COMPREPLY=()
    cur=${COMP_WORDS[COMP_CWORD]}
    prev=${COMP_WORDS[COMP_CWORD-1]}
    cmd=${COMP_WORDS[1]:-}
    if [[ $COMP_CWORD -eq 1 ]]; then
        opts='sketch ani matrix set composite place examples doctor help shuffle dist reverse --help --usage --version --license'
    else
        case "$prev" in
            --metric)
                if [[ $cmd == ani ]]; then
                    opts='best recalibrated ctx-moe ctx-naive p_dist mash aaf mash-if-far aaf-if-far'
                else
                    opts='ctx-moe ctx-naive p_dist mash aaf'
                fi ;;
            --format)
                if [[ $cmd == ani ]]; then opts='detail matrix triangle'
                else opts='full triangle edges clusters dedup-plan'; fi ;;
            --values) opts='distance ani' ;;
            --matrix-format|--keep-matrix-format) opts='tsv phylip' ;;
            --progress) opts='auto on off' ;;
            --dedup-strategy) opts='greedy full-linkage' ;;
            --distance-format) opts='auto ani matrix phylip' ;;
            --rank-by) opts='wls sparse rank_sum' ;;
            *)
                if [[ $cur == -* ]]; then
                    opts=$("${COMP_WORDS[0]}" "$cmd" --help 2>/dev/null |
                        grep -oE -- '--[[:alnum:]][[:alnum:]_-]*' | sort -u)
                else
                    compopt -o filenames 2>/dev/null || true
                    mapfile -t COMPREPLY < <(compgen -f -- "$cur")
                    return 0
                fi ;;
        esac
    fi
    mapfile -t COMPREPLY < <(compgen -W "$opts" -- "$cur")
    return 0
}
complete -F _kssd3a kssd3a
