clc;clear;
X = [12277,12777,13262,13902,14524,15037,15961,16724,17767,19064]; %X为1行N列得向量
% gm11 = GM11_error(X,10,0)      %误差修正模型
% gmV = GMVerhulst_demo(X,10,0)  %背景值优化的GMVerhulst模型
% gm11 = GM11_optimize1(X,10)    %基于背景值优化
% gm11 = GM11_optimize(X,10)     %基于背景值和初始值优化的GM(1,1)模型