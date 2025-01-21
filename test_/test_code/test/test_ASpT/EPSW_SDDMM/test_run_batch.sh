#!/bin/bash

#批量发送可执行
password_c906="rvboards"
password_c910="licheepi"
set source_dir "$PWD"
c906="root@192.168.1.102:/opt/sddmm"
c910="root@192.168.1.112:/mnt/qfc/sddmm"

function send_files(){
  local source_dir=$1
  local remote_dir=$2
  local password=$3

  for file in "$source_dir"*; do
    if [[ -f "$file" && -x "$file" && ! "$file" =~ \.sh$ && ! "$file" =~ \.cc$ && ! "$file" =~ \.png$ && ! "$file" =~ \.py$ ]]; then
      expect -c "
        spawn scp $file $remote_dir
        expect {
          \"*password:\" {
            send $password\r
            exp_continue
          }
        }"

    fi
  done
}

# 调用脚本，并传递根目录、远程目录和密码
send_files "$source_dir" "$c906" "$password_c906"
send_files "$source_dir" "$c910" "$password_c910"