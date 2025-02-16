#This file creates the plots of the paper "Admissible online closed testing must employ e-values"

#Load the simulated data that was generated with the file "generate_data.R"

load("results/TD_data.rda")
load("results/TD_data_taus.rda")

library(ggplot2)

col=c( "aquamarine", "purple", "cornflowerblue",  "darksalmon", "orange", "red", "limegreen", "grey")


#Creates Figure 4 of the paper "Admissible online closed testing must employ e-values"

plot_taus=ggplot(taus_df, aes(x=index, y=taus)) + 
  geom_line(size = 1.2) +
  geom_line(size = 1.2, aes(y=pi_A), col="grey", linetype = "dashed") +
  facet_grid(mu_A ~ pi_A_name)+
  xlab("Hypothesis index")+
  ylab("Estimated proportion of false hypotheses")+
  scale_x_continuous(breaks = seq(10, 90, 20), limits=c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.2), limits=c(0, 1), expand = c(0.02, 0)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y = element_text(size=20), axis.title.x = element_text(size=20),
        legend.text=element_text(size=20), legend.title=element_text(size=20), strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15)) 

ggsave("results/plot_taus.pdf", plot=plot_taus, width=12, height=8)


#Creates Figure 5 of the paper "Admissible online closed testing must employ e-values"

df_KR=TDP_df[TDP_df$Procedure %in% c("online-simple", "closed online-simple", "admissible online-simple", "true proportion"), ]
col_KR=c("aquamarine", "purple", "cornflowerblue", "grey")


line_types = c("online-simple" = "dotdash", 
               "closed online-simple" = "dashed", 
               "admissible online-simple" = "dotted", 
               "true proportion" = "solid")


plot_KR=ggplot(df_KR, aes(x=index, y=TDP, colour = Procedure, linetype = Procedure)) + 
  geom_line(size = 1.2, aes(colour = Procedure)) +
  facet_grid(mu_A ~ pi_A)+
  scale_colour_manual(values = col_KR)+
  scale_linetype_manual(values = line_types)+
  xlab("Hypothesis index")+
  ylab("TDP bound")+
  scale_x_continuous(breaks = seq(100, 900, 200), limits=c(0, 1000), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.2), limits=c(0, 1), expand = c(0.02, 0)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y = element_text(size=20), axis.title.x = element_text(size=20),
        legend.text=element_text(size=20), legend.title=element_text(size=20), strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15)) 

ggsave("results/plot_KR.pdf", plot=plot_KR, width=12, height=8)


#Creates Figure 6 of the paper "Admissible online closed testing must employ e-values"

df_GRO=TDP_df[TDP_df$Procedure %in% c("GRO", "hedged GRO", "boosted GRO", "true proportion"), ]
col_GRO=c("darksalmon", "orange", "red", "grey")

line_types = c("GRO" = "dotdash", 
               "hedged GRO" = "dashed", 
               "boosted GRO" = "dotted", 
               "true proportion" = "solid")

plot_GRO=ggplot(df_GRO, aes(x=index, y=TDP, colour = Procedure, linetype = Procedure)) + 
  geom_line(size = 1.2, aes(colour = Procedure)) +
  facet_grid(mu_A ~ pi_A)+
  scale_colour_manual(values=col_GRO)+
  scale_linetype_manual(values = line_types)+
  xlab("Hypothesis index")+
  ylab("TDP bound")+
  scale_x_continuous(breaks = seq(100, 900, 200), limits=c(0, 1000), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.2), limits=c(0, 1), expand = c(0.02, 0)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y = element_text(size=20), axis.title.x = element_text(size=20),
        legend.text=element_text(size=20), legend.title=element_text(size=20), strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

ggsave("results/plot_GRO.pdf", plot=plot_GRO, width=12, height=8)


#Creates Figure 7 of the paper "Admissible online closed testing must employ e-values"

df_best=TDP_df[TDP_df$Procedure %in% c("admissible online-simple", "boosted GRO", "calibrated", "true proportion"), ]
col_best=c("cornflowerblue", "red", "limegreen", "grey")

line_types = c("calibrated" = "dotdash", 
               "admissible online-simple" = "dashed", 
               "boosted GRO" = "dotted", 
               "true proportion" = "solid")

plot_best=ggplot(df_best, aes(x=index, y=TDP, colour = Procedure, linetype = Procedure)) + 
  geom_line(size = 1.2, aes(colour = Procedure)) +
  facet_grid(mu_A ~ pi_A)+
  scale_colour_manual(values=col_best)+
  scale_linetype_manual(values = line_types)+
  xlab("Hypothesis index")+
  ylab("TDP bound")+
  scale_x_continuous(breaks = seq(100, 900, 200), limits=c(0, 1000), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.2), limits=c(0, 1), expand = c(0.02, 0)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y = element_text(size=20), axis.title.x = element_text(size=20),
        legend.text=element_text(size=20), legend.title=element_text(size=20), strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

ggsave("results/plot_best.pdf", plot=plot_best, width=12, height=8)

#Creates Figure S.1 of the paper "Admissible online closed testing must employ e-values"

df_m_os=TDP_df[TDP_df$Procedure %in% c("m-online-simple", "u-online-simple", "true proportion"), ]
col_m_os=c("grey", "purple", "cornflowerblue")

line_types = c("m-online-simple" = "dashed", 
               "u-online-simple" = "dotted", 
               "true proportion" = "solid")

plot_m_os=ggplot(df_m_os, aes(x=index, y=TDP, colour = Procedure, linetype = Procedure)) + 
  geom_line(size = 1.2, aes(colour = Procedure)) +
  facet_grid(mu_A ~ pi_A)+
  scale_colour_manual(values=col_m_os)+
  scale_linetype_manual(values = line_types)+
  xlab("Hypothesis index")+
  ylab("TDP bound")+
  scale_x_continuous(breaks = seq(100, 900, 200), limits=c(0, 1000), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.2), limits=c(0, 1), expand = c(0.02, 0)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y = element_text(size=20), axis.title.x = element_text(size=20),
        legend.text=element_text(size=20), legend.title=element_text(size=20), strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15)) 

ggsave("results/plot_m-os.pdf", plot=plot_m_os, width=12, height=8)

#Creates Figure S.2 of the paper "Admissible online closed testing must employ e-values"

df_Freed=TDP_df[TDP_df$Procedure %in% c("m-online-Freedman", "u-online-Freedman", "true proportion"), ]
col_Freed=c("grey", "red", "orange")

line_types = c("m-online-Freedman" = "dashed", 
               "u-online-Freedman" = "dotted", 
               "true proportion" = "solid")

plot_Freed=ggplot(df_Freed, aes(x=index, y=TDP, colour = Procedure, linetype = Procedure)) + 
  geom_line(size = 1.2, aes(colour = Procedure)) +
  facet_grid(mu_A ~ pi_A)+
  scale_colour_manual(values=col_Freed)+
  scale_linetype_manual(values = line_types)+
  xlab("Hypothesis index")+
  ylab("TDP bound")+
  scale_x_continuous(breaks = seq(100, 900, 200), limits=c(0, 1000), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.1, 0.9, 0.2), limits=c(0, 1), expand = c(0.02, 0)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y = element_text(size=20), axis.title.x = element_text(size=20),
        legend.text=element_text(size=20), legend.title=element_text(size=20), strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

ggsave("results/plot_m-Freed.pdf", plot=plot_Freed, width=12, height=8)