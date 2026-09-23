#!/bin/sh
# Keep the header timestamp stable unless the recorded build identity changes.
set -eu
root=$1
output=$2
flags=$3
version=$(tr -d '\r\n' < "$root/VERSION")
commit=unknown
state=unknown
snapshot=unknown
checkout=unknown
checkout_state=unknown
if [ -e "$root/.git" ] && git -C "$root" rev-parse --verify HEAD >/dev/null 2>&1; then
    checkout=$(git -C "$root" rev-parse HEAD)
    checkout_state=clean
    if [ -n "$(git -C "$root" status --porcelain --untracked-files=normal)" ]; then checkout_state=dirty; fi
fi
if [ -f "$root/SOURCE_COMMIT" ]; then
    commit=$(awk -F '\t' '$1 == "mother_repo_commit" { print $2 }' "$root/SOURCE_COMMIT")
    state=$(awk -F '\t' '$1 == "mother_repo_status" { print $2 }' "$root/SOURCE_COMMIT")
    snapshot=$(awk -F '\t' '$1 == "source_content_sha256" { print $2 }' "$root/SOURCE_COMMIT")
else
    commit=$checkout
    state=$checkout_state
fi
escape_c() { printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'; }
tmp=$output.tmp.$$
trap 'rm -f "$tmp"' EXIT HUP INT TERM
{
    printf '#define KSSD_VERSION "%s"\n' "$(escape_c "$version")"
    printf '#define KSSD_SOURCE_COMMIT "%s"\n' "$(escape_c "${commit:-unknown}")"
    printf '#define KSSD_SOURCE_STATUS "%s"\n' "$(escape_c "${state:-unknown}")"
    printf '#define KSSD_SOURCE_SNAPSHOT "%s"\n' "$(escape_c "${snapshot:-unknown}")"
    printf '#define KSSD_CHECKOUT_COMMIT "%s"\n' "$(escape_c "$checkout")"
    printf '#define KSSD_CHECKOUT_STATUS "%s"\n' "$(escape_c "$checkout_state")"
    printf '#define KSSD_BUILD_FLAGS "%s"\n' "$(escape_c "$flags")"
} > "$tmp"
if ! cmp -s "$tmp" "$output"; then mv "$tmp" "$output"; fi
