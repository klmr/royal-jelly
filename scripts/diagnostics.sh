_prettify() {
    local id="$(basename "${1%%.*}")"
    id="${id%%Log}"
    local type="$(basename "$(dirname "$1")")"
    echo -n "$type" $'\t' "$id" $'\t'
}

read-length() {
    local file
    for file in data/mapped/*/*final.out
    do
        _prettify "$file"
        awk '/Average input read length/ {print $6/2}' "$file"
    done
}

mapped-fraction() {
    local file
    for file in data/quant/*/*.summary
    do
        _prettify "$file"
        awk -F$'\t' '
            /Unassigned_/ {u=u+$2}
            /Assigned/ {a=$2}
            END {printf("%3.0f%%\n", a / (a + u) * 100)}' \
            "$file"
    done
}

mapped-reads() {
    local file
    for file in data/quant/*/*.summary
    do
        _prettify "$file"
        awk -F$'\t' '
            /Unassigned_/ {u=u+$2}
            /Assigned/ {a=$2}
            END {print a, u}' \
            "$file"
    done
}

insert-size() {
    local file
    for file in data/isize/*/*.tsv
    do
        _prettify "$file"
        grep '^MEDIAN_INSERT_SIZE' --after-context 1 "$file" | tail -1 | cut -f1
    done
}

if [[ $# == 0 ]]; then
    echo 'Available commands:'
    echo
    declare -F | sed 's/.*[[:space:]]\([^[:space:]]*\)$/\1/' | grep -v '^_'
    exit 0
fi

"$@"
