library(deSolve)
library(ggplot2)
library(patchwork)
library(pracma)#找峰值
library(scales)


### 步骤2：定义微分方程
growth_model <- function(t, m, parameters) {
  with(as.list(parameters), {
    dm_dt <- a * m^(alpha) - (b+p/max_m) * m #west 方程
    #dm_dt <- g * m * (1-m/max_m)*(1-p/p_max)#deng方程
    return(list(dm_dt))
  })
}
# 初始条件
set.seed(123)  # 设置随机数种子
random_numbers1 <- rnorm(100, mean = 50, sd = 10)# 初始生物量
random_numbers2 <- rnorm(100, mean = 50, sd = 5)# 初始生物量
random_numbers3 <- rnorm(100, mean = 50, sd = 2)# 初始生物量

random_numbers1=sort(random_numbers1)#升序排序
random_numbers2=sort(random_numbers2)#升序排序
random_numbers3=sort(random_numbers3)#升序排序
rand_tot=data.frame(tot_sd10=random_numbers1,tot_sd5=random_numbers2,tot_sd2=random_numbers3)
totmass_v=c(random_numbers1,random_numbers2,random_numbers3)

# 参数
g=0.1
a <- 0.1  # 用于生长的能量分配系数
b <- 0.01  # 用于维持的能量分配系数
###################
#without competition west's equation
p=0
p_max=10
alpha=0.75
theta=0.5#临界阈值
beta=1
max_m=1000000
leafratio=matrix(nrow=100,ncol=100);int_lf=end_lf=matrix(nrow=100,ncol=3)
#df=list()
for (j in 1:3){
for (i in 1:100){
  m0=rand_tot[i,j]
# 时间序列
times <- seq(1, 100, by = 1)  # 从1到100，步长为1
### 步骤4：求解微分方程
out <- ode(y = m0, times = times, func = growth_model, parms = c(a = a, b = b,p=p))
total=out[,2]
  if(i<theta*100){leafratio[,i]=beta*(1+p/p_max)*total^(alpha-1)}#随竞争增加叶比重增加
  if(i>=theta*100){leafratio[,i]=beta*(1-p/p_max)*total^(alpha-1)}#随竞争增加叶比重减小
}
  
  int_lf[,j]=round(leafratio[10,],4)#t=10时刻
  end_lf[,j]=round(leafratio[100,],4)#最终时刻 放附件
}
int_df=c(int_lf[,1],int_lf[,2],int_lf[,3])
end_df=c(end_lf[,1],end_lf[,2],end_lf[,3])
class=c(rep("sd = 10",100),rep("sd = 5",100),rep("sd = 2",100))
df=data.frame(int_df=int_df,end_df=end_df,class=class,totmass_v=totmass_v)

totalmass=ggplot(df,aes(x = totmass_v,fill = class)) + geom_density(alpha=0.5)+
  labs(title = "Initial state", x = "Total biomass",y="Probability denstiy") +
  annotate("text",x=20, y=0.29, label="A",angle = 0,size=15,family="sans")+
  theme(axis.text=element_text(size=25,family="serif",colour = "black"),
        axis.text.y=element_text(size=25,family="serif",angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
        axis.text.x=element_text(size=25,family="serif",angle=0,colour = "black",vjust = 0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.85,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 23,color='black'),
        legend.title = element_text(size = 16,color=NA),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2.5,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm"))+
  scale_x_continuous(breaks=seq(20,80,by=20),limits=c(20,80))+
  scale_y_continuous(breaks=seq(0,0.3,by=0.1),limits=c(0,0.3))


nocom_t10=ggplot(df,aes(x = int_df,fill = class)) + geom_density(alpha=0.5)+
  labs(title = "p = 0; t = 10", x = "Leaf mass ratio",y="Probability denstiy") +
  annotate("text",x=0.25, y=140, label="B",angle = 0,size=15,family="sans")+
  theme(axis.text=element_text(size=25,family="serif",colour = "black"),
        axis.text.y=element_text(size=25,family="serif",angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
        axis.text.x=element_text(size=25,family="serif",angle=0,colour = "black",vjust = 0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.85,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 23,color='black'),
        legend.title = element_text(size = 16,color=NA),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2.5,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm"))+
  scale_x_continuous(breaks=seq(0.25,0.4,by=0.05),limits=c(0.25,0.4))+
  scale_y_continuous(breaks=seq(0,150,by=50),limits=c(0,150))

nocom_end=ggplot(df,aes(x = end_df,fill = class)) + geom_density(alpha=0.5)+
  labs(title = "p = 0; t = 100", x = "Leaf mass ratio",y="Probability denstiy") +
  annotate("text",x=0.21, y=485, label="A",angle = 0,size=15,family="sans")+
  theme(axis.text=element_text(size=25,family="serif",colour = "black"),
        axis.text.y=element_text(size=25,family="serif",angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
        axis.text.x=element_text(size=25,family="serif",angle=0,colour = "black",vjust = 0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.85,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 23,color='black'),
        legend.title = element_text(size = 16,color=NA),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2.5,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm"))+
  scale_x_continuous(breaks=seq(0.21,0.27,by=0.02),limits=c(0.21,0.27))+
  scale_y_continuous(breaks=seq(0,500,by=500/4),limits=c(0,500))
##########################
#with competition west's equation
p=2
leafratio=matrix(nrow=100,ncol=100);int_lf=end_lf=matrix(nrow=100,ncol=3)
for (j in 1:3){
  for (i in 1:100){
    m0=rand_tot[i,j]
    # 时间序列
    times <- seq(1, 100, by = 1)  # 从1到100，步长为1
    ### 步骤4：求解微分方程
    out <- ode(y = m0, times = times, func = growth_model, parms = c(a = a, b = b,p=p))
    total=out[,2]
    if(i<theta*100){leafratio[,i]=beta*(1+p/p_max)*total^(alpha-1)}#随竞争增加叶比重增加
    if(i>=theta*100){leafratio[,i]=beta*(1-p/p_max)*total^(alpha-1)}#随竞争增加叶比重减小
  }
  
  int_lf[,j]=round(leafratio[10,],4)#t=10时刻
  end_lf[,j]=round(leafratio[100,],4)#最终时刻 放附件
}
int_df=c(int_lf[,1],int_lf[,2],int_lf[,3])
end_df=c(end_lf[,1],end_lf[,2],end_lf[,3])
class=c(rep("sd = 10",100),rep("sd = 5",100),rep("sd = 2",100))
df=data.frame(int_df=int_df,end_df=end_df,class=class)

com_t10=ggplot(df,aes(x = int_df,fill = class)) + geom_density(alpha=0.5)+
  labs(title = "p = 2; t = 10", x = "Leaf mass ratio",y="Probability denstiy") +
  annotate("text",x=0.1, y=14.5, label="C",angle = 0,size=15,family="sans")+
  theme(axis.text=element_text(size=25,family="serif",colour = "black"),
        axis.text.y=element_text(size=25,family="serif",angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
        axis.text.x=element_text(size=25,family="serif",angle=0,colour = "black",vjust = 0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.85,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 23,color='black'),
        legend.title = element_text(size = 16,color=NA),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2.5,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm"))+
  scale_x_continuous(breaks=seq(0.1,0.7,by=0.2),limits=c(0.1,0.7))+
  scale_y_continuous(breaks=seq(0,15,by=5),limits=c(0,15))

com_end=ggplot(df,aes(x = end_df,fill = class)) + geom_density(alpha=0.5)+
  labs(title = "p = 2; t = 100", x = "Leaf mass ratio",y="Probability denstiy") +
  annotate("text",x=0.1, y=19, label="B",angle = 0,size=15,family="sans")+
  theme(axis.text=element_text(size=25,family="serif",colour = "black"),
        axis.text.y=element_text(size=25,family="serif",angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
        axis.text.x=element_text(size=25,family="serif",angle=0,colour = "black",vjust = 0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.8,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 23,color='black'),
        legend.title = element_text(size = 16,color=NA),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2.5,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm"))+
  scale_x_continuous(breaks=seq(0.1,0.4,by=0.1),limits=c(0.1,0.4))+
  scale_y_continuous(breaks=seq(0,20,by=5),limits=c(0,20))

#####################################
#sensitive analyses in parameter space
alphav=c(0.6,0.7,0.8,0.9,1) #scaling exponent
pv=seq(0,4,1)  #competition
scaling_exponent=competition=total_distance=average_distance=NA
kk=1
for (ij1 in 1:length(alphav)){
  for (ij2 in 1:length(pv)){
    alpha=alphav[ij1]
    p=pv[ij2]
    
    leafratio=matrix(nrow=100,ncol=100)
    for (i in 1:100){
      m0=random_numbers1[i]
      # 时间序列
      times <- seq(1, 100, by = 1)  # 从1到100，步长为1
      ### 步骤4：求解微分方程
      out <- ode(y = m0, times = times, func = growth_model, parms = c(a = a, b = b,p=p))
      total=out[,2]
      if(i<theta*100){leafratio[,i]=beta*(1+p/p_max)*total^(alpha-1)}#随竞争增加叶比重增加
      if(i>=theta*100){leafratio[,i]=beta*(1-p/p_max)*total^(alpha-1)}#随竞争增加叶比重减小
    }
    
    lf_t10=round(leafratio[10,],4)#t=10时刻
    pdf_obj=density(lf_t10)
    peak_indices=findpeaks(pdf_obj$y, npeaks=3,minpeakheight  = 0.001, minpeakdistance = 1,sortstr=TRUE)
    peak_values <- pdf_obj$x[peak_indices[,2]]
    peak_densities <- peak_indices[,1]
    peak_df=data.frame(peak_values=peak_values,peak_densities=peak_densities,
                       scaling_exponent=rep(alpha,length(peak_densities)),competition=rep(p,length(peak_densities)))
    #peakT=rbind(peakT,peak_df)
    # 计算每个值之间的距离
    distance_matrix <- outer(peak_values, peak_values, function(ii, jj) abs(ii - jj))
    # 计算所有距离的总和
    total_distance[kk] <- sum(distance_matrix)
    # 计算值对的总数（不包括自身）
    num_pairs = length(peak_values)
    if(num_pairs>1){
      total_pairs = num_pairs * (num_pairs - 1)
      # 计算平均距离
      average_distance[kk]=sum(distance_matrix)/total_pairs
    }else{
      average_distance[kk]=sum(distance_matrix)
    }
    scaling_exponent[kk]=alpha
    competition[kk]=p
    kk=kk+1
  }
}

distance_df=data.frame(total_distance=na.omit(total_distance),average_distance=na.omit(average_distance),
                       scaling_exponent=na.omit(scaling_exponent),competition=na.omit(competition))
sense_totdis=ggplot(distance_df,aes(x=scaling_exponent,y=competition,fill=total_distance))+geom_tile(position = "identity")+
  annotate("text",x=0.55, y=4.3, label="C",angle = 0,size=15,family="sans")+
  scale_fill_gradientn(colors = c("#08306B",  "#9ECAE1", "#C6DBEF", "#DEEBF7"), 
                       values = scales::rescale(c(0.0,0.1,0.2,0.3)), 
                       breaks = c(0.0, 0.8, 1.59), 
                       labels = c("0.0", "0.8", "1.6"))+
labs(title = "Total distance", x = "Scaling exponent", y = "Light stress intensity (p)", fill = "") +
  theme(axis.text=element_text(size=25,family="serif",colour = "black"),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom", 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 20,color='black'),
        legend.title = element_text(size = 25,color='black'),
        legend.direction="horizontal", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm"))

sense_avgdis=ggplot(distance_df,aes(x=scaling_exponent,y=competition,fill=average_distance))+geom_tile(position = "identity")+
  scale_fill_gradientn(colors = c("#08306B",  "#9ECAE1", "#C6DBEF", "#DEEBF7"), 
                       values = scales::rescale(c(0.05, 0.1, 0.2, 0.31)), 
                       breaks = c(0.0, 0.4,0.797), 
                       labels = c("0.0", "0.4","0.8"))+
annotate("text",x=0.55, y=4.3, label="D",angle = 0,size=15,family="sans")+
  labs(title = "Average distance", x = "Scaling exponent", y = "Light stress intensity (p)", fill = "") +
  theme(axis.text=element_text(size=25,family="serif",colour = "black"),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.5, colour = "black", linetype=1),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 1, color="black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom", 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 20,color='black'),
        legend.title = element_text(size = 25,color='black'),
        legend.direction="horizontal", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm"))

#table(round(leafratio[1,],4))
#table(round(leafratio[100,],4))

#ggsave("West model_fig1.pdf",(totalmass+nocom_t10)/(com_t10+ sense_avgdis),width = 40, height = 40, units = "cm", dpi = 300) 
#ggsave("West model_fig2.pdf",(nocom_end+ com_end+sense_totdis),width = 60, height = 20, units = "cm", dpi = 300) 
