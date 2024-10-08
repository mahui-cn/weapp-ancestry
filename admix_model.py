# -*- coding: utf-8 -*-
import os
import csv
import numpy as np
import copy


# 祖源模型类
class AdmixModel:

    # 祖源模型数据的相对路径
    __admixModelPath = ""

    # 祖源模型信息
    __admixModelInfo = [
        {
            "key": "wbbc",
            "name": "Westlake BioBank for Chinese K4",
            "name_cn": "西湖中国 K4",
            "desc": "Westlake BioBank for Chinese is a prospective study based on the Chinese population, with the main purpose of better understanding the impact of genetic and environmental factors on growth and development from young to old age. In this study, 10,376 men and women were recruited from 29 provinces in China (excluding Beijing, Shanghai, Xinjiang, Tibet, Hong Kong, Macao and Taiwan) for Whole Genome Sequencing, with an average sequencing depth of 13.9x, a total of 81.5 million SNPs were sequenced, accounting for 99.77% of the human genome, and 4480 samples were selected for subsequent analysis after quality control. PCR was performed based on the SNP mutation rate, and all samples were divided into four Han subgroups (north, central, south and Lingnan), corresponding to the geographical boundaries of the Qinling-Huaihe and Nanling mountain, respectively. This admixture model is available for East, North and Southeast Asian only.",
            "desc_cn": "西湖（中国）生物库项目是一项基于中国人群的前瞻性研究，其主要目的是更好地了解遗传和环境因素对从年轻到老年的生长发育的影响。本研究从中国29个省（不含北京、上海、新疆、西藏和港澳台地区）广泛招募了10,376个男性和女性进行WGS全基因组测序，平均测序深度为13.9x，共计测序8150万个SNP，占比人类基因组99.77%，经过质控筛选出4480个样本进行后续分析。根据SNP突变率进行PCA主成分分析，以及秦岭-淮河和南岭山脉的地理边界，对所有样本划分为四个汉族子群（北部、中部、南部和岭南）。此祖源模型只适用于东亚、北亚、东南亚人群。",
            "admix": [
                {
                    "name": "North Han Chinese",
                    "name_cn": "中国北方汉族",
                    "desc": "North Han Chinese include people in Shandong, Shanxi, Tianjin, Hebei, Henan, Inner Mongolia, Heilongjiang, Jilin, Liaoning, Shaanxi, Gansu, Ningxia and Qinghai.",
                    "desc_cn": "中国北方汉族包括山东、山西、天津、河北、河南、内蒙古、黑龙江、吉林、辽宁、陕西、甘肃、宁夏和青海。",
                },
                {
                    "name": "Central Han Chinese",
                    "name_cn": "华中汉族",
                    "desc": "Central Han Chinese include people in Anhui and Jiangsu.",
                    "desc_cn": "华中汉族包括安徽和江苏。",
                },
                {
                    "name": "South Han Chinese",
                    "name_cn": "中国南方汉族",
                    "desc": "South Han Chinese include people in Chongqing, Fujian, Guizhou, Hubei, Hunan, Jiangxi, Sichuan, Yunnan, and Zhejiang.",
                    "desc_cn": "中国南方汉族包括重庆、福建、贵州、湖北、湖南、江西、四川、云南和浙江。",
                },
                {
                    "name": "Lingnan Han Chinese",
                    "name_cn": "岭南汉族",
                    "desc": "Lingnan Han Chinese include people in Guangxi, Guangdong and Hainan",
                    "desc_cn": "岭南汉族包括广西、广东和海南。",
                },
            ],
        }
    ]

    @property
    def AdmixModelPath(self):
        return self.__admixModelPath

    @property
    def AdmixModelInfo(self):
        return self.__admixModelInfo

    def __init__(self, admixModelPath="model"):
        self.__admixModelPath = admixModelPath

    # 根据祖源模型名字获取模型信息，深拷贝不修改源对象
    def get_model_info(self, model_name):
        for model in self.__admixModelInfo:
            if model["key"].lower() == model_name.lower():
                return copy.deepcopy(model)
        raise Exception("不支持祖源模型：" + model_name)

    # 获取祖源模型的SNP和Frequency数据
    def get_model_data(self, model_name):
        # 祖源模型数据源文件
        snp_file_name = os.path.join(self.__admixModelPath, model_name + ".alleles")
        frequency_file_name = os.path.join(self.__admixModelPath, model_name + ".F")

        if not os.access(snp_file_name, os.F_OK):
            raise FileNotFoundError("祖源模型SNP文件不存在：" + snp_file_name)

        if not os.access(frequency_file_name, os.F_OK):
            raise FileNotFoundError("祖源模型SNP频率文件不存在：" + frequency_file_name)

        # 读取SNP的rsid, minor, major突变
        rsid = []
        minor_geno = []
        major_geno = []
        with open(
            snp_file_name,
            "r",
        ) as snp_file:
            for row in csv.reader(snp_file, delimiter=" "):
                if len(row) == 3:
                    rsid.append(row[0])
                    minor_geno.append(row[1])
                    major_geno.append(row[2])

        # 读取SNP major突变频率
        frequency = []
        with open(
            frequency_file_name,
            "r",
        ) as frequency_file:
            for row in csv.reader(frequency_file, delimiter=" "):
                frequency.append([float(f) for f in row])

        if len(rsid) == 0:
            raise Exception("祖源模型“{}”数据文件中没有数据".format(model_name))

        if not len(rsid) == len(minor_geno) == len(major_geno) == len(frequency):
            raise Exception(
                "祖源模型“{}”的alleles和frequency数据行数不一致".format(model_name)
            )

        return (
            np.array(rsid),
            np.array(major_geno),
            np.array(minor_geno),
            np.array(frequency),
        )
