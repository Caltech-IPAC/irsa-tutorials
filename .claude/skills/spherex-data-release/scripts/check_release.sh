#!/usr/bin/env bash
# Report what is live for a SPHEREx release on IRSA's pre-release and production hosts.
# Usage: check_release.sh <release>   e.g. check_release.sh dr2
# The irsadev/irsatest checks need the IRSA network or VPN.
set -u
REL="${1:?usage: check_release.sh <release, e.g. dr2>}"

sia_collections() {
    local q="SELECT DISTINCT obs_collection FROM ivoa.obscore WHERE obs_collection LIKE 'spherex%'"
    curl -s -m 60 --data-urlencode "QUERY=$q" -d FORMAT=csv -d LANG=ADQL "$1/TAP/sync" \
        | tail -n +2 | tr '\n' ' '
}

tap_rows() {
    local q="SELECT TOP 1 a.uri FROM spherex.artifact a WHERE a.uri LIKE 'ibe/data/spherex/$REL/%'"
    curl -s -m 120 --data-urlencode "QUERY=$q" -d FORMAT=csv -d LANG=ADQL "$1/TAP/sync" | sed -n 2p
}

# Ask SIA, not TAP: ivoa.obscore leaves cloud_access blank even where SIA (what the notebooks
# use) fills it. Search at the position of one of the release's own images. IRSA's SIA CSV
# doesn't quote the JSON, so grep it out instead of parsing columns.
cloud_access() {
    local q="SELECT TOP 1 s_ra, s_dec FROM ivoa.obscore WHERE obs_collection = 'spherex_$REL'"
    local pos ra dec row
    pos=$(curl -s -m 60 --data-urlencode "QUERY=$q" -d FORMAT=csv -d LANG=ADQL "$1/TAP/sync" | sed -n 2p)
    [ -z "$pos" ] && { echo "(no spherex_$REL rows)"; return; }
    IFS=, read -r ra dec <<<"$pos"
    row=$(curl -s -m 120 "$1/SIA?COLLECTION=spherex_$REL&POS=CIRCLE+$ra+$dec+0.0003&MAXREC=1&RESPONSEFORMAT=csv" | sed -n 2p)
    [ -z "$row" ] && { echo "(no SIA rows)"; return; }
    grep -o '{"aws": {[^}]*}}' <<<"$row" || echo "(blank)"
}

# Subfolders under <release>/ in a bucket ("(empty)" if none, or the S3 error code).
s3_prefixes() {
    local xml out
    xml=$(curl -s -m 60 "https://$1.s3.us-east-1.amazonaws.com/?list-type=2&delimiter=/&prefix=$REL/")
    if grep -q '<Error>' <<<"$xml"; then
        grep -o '<Code>[^<]*</Code>' <<<"$xml" | sed 's/<[^>]*>//g'
        return
    fi
    out=$(grep -o '<CommonPrefixes><Prefix>[^<]*</Prefix>' <<<"$xml" | sed 's/<[^>]*>//g' | tr '\n' ' ')
    echo "${out:-(empty)}"
}

http_code() { curl -s -m 30 -o /dev/null -I -w '%{http_code}' "$1"; }

for host in https://irsadev.ipac.caltech.edu https://irsatest.ipac.caltech.edu https://irsa.ipac.caltech.edu; do
    echo "== $host"
    echo "  SIA collections:   $(sia_collections "$host")"
    uri=$(tap_rows "$host")
    echo "  $REL TAP row:       ${uri:-(none)}"
    echo "  cloud_access (spherex_$REL): $(cloud_access "$host")"
    # IBE has no directory listing, so check that the TAP row's file is served.
    [ -n "$uri" ] && echo "  IBE file for that row: HTTP $(http_code "$host/$uri")"
done

echo "== S3"
echo "  staging nasa-irsa-spherex-$REL-staging: $(s3_prefixes "nasa-irsa-spherex-$REL-staging")"
echo "  public  nasa-irsa-spherex:              $(s3_prefixes nasa-irsa-spherex)"

echo "== Docs"
for f in SPHEREx_Expsupp.pdf SPHEREx_Expsupp_QR.pdf; do
    echo "  $f: $(http_code "https://irsa.ipac.caltech.edu/data/SPHEREx/docs/$f")"
done
