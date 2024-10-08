# -*- coding: utf-8 -*-
import numpy as np
import scipy.optimize as optimize
from admix_model import *


# 祖源计算器类
class AdmixCalculator:

    # 计算用户每个alleles等位基因相对于祖源模型的major和minor基因型的突变总次数
    def __get_geno_stat(self, user_genome, rsid, major_geno, minor_geno):

        major_geno_count = []
        minor_geno_count = []

        # allele点位匹配计数
        match_count = 0

        for i in range(len(rsid)):
            major_geno_count.append(0)
            minor_geno_count.append(0)
            if rsid[i] in user_genome:
                match_count += 1
                user_geno = user_genome.get(rsid[i])
                if len(user_geno) == 2:
                    allele1 = user_geno[0]
                    allele2 = user_geno[1]
                elif len(user_geno) == 1:
                    allele1 = user_geno[0]
                    allele2 = user_geno[0]
                else:
                    allele1 = "-"
                    allele2 = "-"

                if allele1 == major_geno[i]:
                    major_geno_count[i] += 1

                if allele2 == major_geno[i]:
                    major_geno_count[i] += 1

                if allele1 == minor_geno[i]:
                    minor_geno_count[i] += 1

                if allele2 == minor_geno[i]:
                    minor_geno_count[i] += 1

        # 匹配的allele点位在祖源模型中的占比
        model_match_ratio = match_count / len(rsid) if len(rsid) > 0 else 0

        # 匹配的allele点位在用户数据中的占比
        user_match_ratio = match_count / len(user_genome) if len(user_genome) > 0 else 0

        return (
            np.array(major_geno_count),
            np.array(minor_geno_count),
            model_match_ratio,
            user_match_ratio,
        )

    # 损失函数
    def __loss_func(self, major_geno_count, minor_geno_count, frequency, admix_ratio):
        # 计算用户所有等位基因发生major和minor突变频率的加权平均值
        major_frq_mean = np.dot(frequency, admix_ratio)
        minor_frq_mean = np.dot(1 - frequency, admix_ratio)

        # 计算用户所有等位基因发生major和minor突变的联合概率
        major_pr = np.dot(
            major_geno_count, np.log(np.where(major_frq_mean > 0, major_frq_mean, 1))
        )
        minor_pr = np.dot(
            minor_geno_count, np.log(np.where(minor_frq_mean > 0, minor_frq_mean, 1))
        )

        # 返回负值作为损失值
        return -(major_pr + minor_pr)

    # 根据用户基因数据，祖源模型，使用极大似然估计计算用户的祖源成分
    def calc_admix(self, user_genome, admix_model_name, opt_tol=1e-4):
        try:
            # 获取指定的祖源模型信息
            admixModel = AdmixModel()
            admix_info = admixModel.get_model_info(admix_model_name)
            admix_count = len(admix_info["admix"])
            if admix_count < 1:
                raise Exception("祖源模型“{}”没有人群成分".format(admix_model_name))

            # 获取指定祖源模型的rsid, alleles数据源
            rsid, major_geno, minor_geno, frequency = admixModel.get_model_data(
                admix_model_name
            )

            # 祖源模型的SNP数量
            snp_count = len(rsid)

            if admix_count != len(frequency[0]):
                raise Exception(
                    "祖源模型“{}”的信息文件和数据文件的人群成分数量不一致".format(
                        admix_model_name
                    )
                )

            # 计算用户基因相对于祖源模型的突变情况
            major_geno_count, minor_geno_count, user_match_ratio, model_match_ratio = (
                self.__get_geno_stat(user_genome, rsid, major_geno, minor_geno)
            )

            # 损失函数的参数初始值
            initial_guess = np.ones(admix_count) / admix_count

            # 损失函数的参数约束条件
            bounds = tuple((0, 1) for i in range(admix_count))
            constraints = {"type": "eq", "fun": lambda ratio: np.sum(ratio) - 1}

            # 计算用户祖源成分的最优值
            res = optimize.minimize(
                fun=lambda ratio: self.__loss_func(
                    major_geno_count, minor_geno_count, frequency, ratio
                ),
                x0=initial_guess,
                bounds=bounds,
                constraints=constraints,
                tol=opt_tol,
            )

            if res.success:
                admix_ratio = res.x
            else:
                admix_ratio = np.zeros(admix_count)

            # 用户祖源成分值写入祖源模型
            for i, admix in enumerate(admix_info["admix"]):
                admix["ratio"] = admix_ratio[i]

            admix_info["snp_count"] = snp_count
            admix_info["user_match_ratio"] = user_match_ratio
            admix_info["model_match_ratio"] = model_match_ratio

            # 对祖源成分比例降序排序
            admix_info["admix"].sort(key=lambda admix: admix["ratio"], reverse=True)

            return admix_info
        except Exception as e:
            raise e
