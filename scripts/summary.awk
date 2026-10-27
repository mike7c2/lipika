/summary: \{/ { on = 1; next }
on && /^    \},/ { exit }
on && /^      [a-z_]+: [0-9.]+,$/ { gsub(/[ ,]/, ""); split($0, kv, ":"); printf "  %-16s %s\n", kv[1], kv[2]; next }
on && /^      [a-z_]+: \{$/ { key = $1; sub(/:$/, "", key); getline; gsub(/[ ,]/, ""); sub(/cnt:/, ""); printf "  %-16s %s\n", key, $0 }
