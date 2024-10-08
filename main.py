# -*- coding: utf-8 -*-
import sys
import json
import warnings

# wegene_utils 库会包含在每个应用的环境当中，无需自行打包上传
# 这里提供源代码以供应用完整运行
import wegene_utils

from admix_calculator import *

# import echarts
import matplot

"""
当输入是部分位点时, 基因位点数据以 json 形式输入:
    {"inputs": {"RS671": "AA", "RS12203592": "CA", "format": "wegene_affy_2"}}
当输入全部位点时，全部位点对应的字符串序列会被 gzip 压缩并以 base64 转码，
转码前的数据在 data 域中
    {"inputs": {"data": "xfgakljdflkja...", format: "wegene_affy_2"}}
你需要解码并解压该数据，解压后的字符串如下:
    AACCTACCCCCC...
进一步，你需要利用相应格式的索引文件对每个位点进行解析
我们提供了整个流程的示例代码以及索引文件
"""

warnings.filterwarnings("ignore")

# 从 stdin 读取输入数据
body = sys.stdin.read()

try:
    # 如果输入的数据是全部位点数据，可以使用 wegene_utils的方法进行解析后使用
    inputs = json.loads(body)["inputs"]
    # 使用 wegene_utils 解析以后，数据会被解析成下面这样的 json 格式:
    #   {'rs123': {'genotype': 'AA', 'chromosome': '1', position: '1236'}, ...}
    user_genome = wegene_utils.process_raw_genome_data(inputs)

    # 筛选出用户的常染SNP
    chr_set = {str(chr) for chr in range(1, 23)}
    user_snp_dict = {}
    for rsid in user_genome.keys():
        if (
            "genotype" in user_genome[rsid]
            and "chromosome" in user_genome[rsid]
            and user_genome[rsid]["chromosome"] in chr_set
            and len(user_genome[rsid]["genotype"]) == 2
            and user_genome[rsid]["genotype"][0] in {"A", "T", "G", "C"}
        ):
            user_snp_dict[rsid] = user_genome[rsid]["genotype"]

    if len(user_snp_dict) == 0:
        raise Exception("用户的基因数据没有有效的SNP")

    # 计算祖源成分
    admixCalc = AdmixCalculator()
    admix_info = admixCalc.calc_admix(user_snp_dict, "wbbc")

    # 根据祖源成分生成图
    pie_x = []
    pie_label = []
    for admix in admix_info["admix"]:
        ratio = round(admix["ratio"] * 100, 2)
        if ratio > 0:
            pie_x.append(ratio)
            pie_label.append(admix["name"])

    pie_base64 = matplot.make_pie(pie_x, pie_label) if len(pie_x) > 0 else None

    # 输出HTML
    result = []
    if pie_base64 != None:
        result.append(
            "<img src='data:image/{};base64,{}' class='img-thumbnail mx-auto d-block'/>".format(
                "png", pie_base64
            )
        )
    result.append("<div class='table-responsive'>")
    result.append("<table class='table'>")
    result.append("<thead><tr><th>祖源成分</th><th>比例</th></tr></thead>")
    result.append("<tbody>")
    for admix in admix_info["admix"]:
        result.append(
            "<tr><td>{}</td><td>{}%</td></tr>".format(
                admix["name_cn"], round(admix["ratio"] * 100, 2)
            )
        )
    result.append("</tbody>")
    result.append("</table>")
    result.append("</div>")
    result.append("<ul>")
    for admix in admix_info["admix"]:
        result.append("<li>{}</li>".format(admix["desc_cn"]))
    result.append("</ul>")
    result.append(
        "<hr/><div class='alert alert-info' style='margin-top: 10px;' role='alert'>基于西湖样本库，此祖源计算器从全球基因公司的基因芯片数据中精选{:,}个SNP作为参考数据集，此次计算匹配了您的{:,}个常染色体SNP的{:.2%}。非WeGene用户，可在<a href='http://geneu.xyz' target='_blank'>基因助手GeneU.xyz</a>上传样本使用此祖源计算器。</div>".format(
            admix_info["snp_count"],
            len(user_snp_dict),
            admix_info["user_match_ratio"],
        )
    )

    # 输出给用户的结果只需要通过 print 输出即可，print只可调用一次
    print("".join(result))

except Exception as e:
    # 错误信息需要被从 stderr 中输出，否则会作为正常结果输出
    for msg in e.args:
        sys.stderr.write(msg)
    exit(2)
