#!/bin/bash

# 检查是否安装了fprettify
if ! command -v fprettify
then
    has_fprettify=$(pip list | grep fprettify)
    # 存在fprettify但没有相应的环境变量，导出环境变量
    if [ -z "$has_fprettify" ]; then
        echo "fprettify未安装,尝试自动安装..."
        pip install fprettify
        has_fprettify=$(pip list | grep fprettify)
        if [ -n "$has_fprettify" ]; then
            echo "fprettify安装成功。"
        else
            echo "fprettify安装失败。"
            exit 1
        fi
    fi
    export fprettify_path="$( dirname $(dirname $(dirname $(pip show fprettify | grep "Location:" | awk '{print $2}'))))/bin"
    export PATH=$PATH:$fprettify_path
    # 添加环境变量到~/.bashrc
    if [ -n "$fprettify_path" ]; then
        echo "将fprettify路径添加到bashrc环境变量。"
        echo 'export PATH=$PATH:'"$fprettify_path" >> ~/.bashrc
        echo 'export fprettify=$fprettify_path/fprettify' >> ~/.bashrc
        echo "fprettify环境变量添加成功。"
    fi
    if ! command -v fprettify 
    then
        echo "fprettify环境变量添加失败。"
        exit 1
    fi
fi
echo "fprettify已安装。"

# 获取修改过的文件列表
changed_files=$(git diff --name-only HEAD)

# 对列表中后缀为F90或f90的文件进行格式化
for file in $changed_files; do
    if [[ -f $file && ($file == *.F90 || $file == *.f90) ]]; then
        fprettify -i 2 -w 3 -l 600 --case 2 2 2 2 -e '*.for' -e '*.f' $file
        echo "Formatting $file done!"
    fi
    
    # -i 2 表示缩进为2个空格
    # -w 3 表示自动预设 运算符、打印/读取、加/减/乘/除运算符的空格
    # -l 600 表示每行字符数不超过600
    # --case 2 表示关键字全部大写
    # -e '*.for' -e '*.f' 表示排除.for和.f文件
done