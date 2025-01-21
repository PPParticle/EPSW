#!/bin/bash

#批量发送可执行
set source_dir "$PWD"
c908="root@192.168.1.105:/tmp/data/spmspm"
password_c908=""

function send_files(){
  local source_dir=$1
  local remote_dir=$2
  local password=$3

  for file in "$source_dir"*_xt*; do
    if [[ -f "$file" && -x "$file" && ! "$file" =~ \.sh$ && ! "$file" =~ \.cc$ ]]; then
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
send_files "$source_dir" "$c908" "$password_c908"
