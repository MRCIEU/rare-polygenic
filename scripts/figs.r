library(dplyr)
library(ggplot2)
library(here)

load(file=here("results/run_sim.RData"))

l1 %>% 
    filter(out == "prs", family_rank == 1) %>%
    group_by(gen, what, h2, rho) %>%
    summarise(r=median(r^2)) %>%
    mutate(gen=gen+1, what = case_when(what == "rare" ~ "Null rare variants", what == "common" ~ "Null common variants", what == "prs" ~ "Causal variants")) %>%
    mutate(rholab = paste("rho ==", rho), h2lab = paste("h^2 ==", h2)) %>%
    ggplot(., aes(x=as.factor(gen), y=r, colour=as.factor(what))) + geom_line(aes(group=what)) +
    facet_grid(h2lab ~ rholab, labeller = label_parsed) +
    labs(x="Generation", y=expression("Median" ~ R^2 ~ "between variant and PRS"), colour="Variant") +
    scale_colour_brewer(type="qual") +
    theme_minimal()
ggsave(here("results/run_sim_prs.pdf"), width=7, height=7)

l1 %>% 
    filter(out == "y", family_rank == 1) %>%
    group_by(gen, what, h2, rho) %>%
    summarise(r=median(r^2)) %>%
    mutate(gen=gen+1, what = case_when(what == "rare" ~ "Null rare variants", what == "common" ~ "Null common variants", what == "prs" ~ "Causal variants")) %>%
    mutate(rholab = paste("rho ==", rho), h2lab = paste("h^2 ==", h2)) %>%
    ggplot(., aes(x=as.factor(gen), y=r, colour=as.factor(what))) + geom_line(aes(group=what)) +
    facet_grid(h2lab ~ rholab, labeller = label_parsed) +
    labs(x="Generation", y=expression("Median" ~ R^2 ~ "between variant and trait"), colour="Variant") +
    scale_colour_brewer(type="qual") +
    theme_minimal()
ggsave(here("results/run_sim_trait.pdf"), width=7, height=7)


load(file=here("results/run_sim2.RData"))
temp <- l2 %>% 
    filter(out == "y") %>%
    group_by(family_rank, what, h2, rho, gen) %>%
    summarise(r=median(r^2)) %>%
    mutate(
        rholab = paste("rho ==", rho), h2lab = paste("h^2 ==", h2)
    )

temp <- mutate(temp, 
    what = case_when(
        what == "rare" ~ "Null rare variants", 
        what == "common" ~ "Null common variants", 
        what == "prs" ~ "Causal variants"
    )
)

p1 <- temp %>%
    filter(family_rank != 0.5) %>%
    ggplot(., aes(x=gen, y=r)) + 
    geom_point(aes(colour=as.factor(family_rank))) +
    geom_smooth(aes(colour=as.factor(family_rank)), se=FALSE) +
    facet_grid(. ~ what) +
    labs(x="Generation", y=expression("Median" ~ R^2 ~ "between variant and trait"), linetype="Variant", colour="Phenotypic rank") +
    scale_colour_brewer(type="seq") +
    theme_minimal() +
    theme(axis.text.x= element_text(angle=90, hjust=0.5))
ggsave(p1, file="results/run_sim_rank.pdf", width=14)


load(file=here("results/run_sim3.RData"))
temp <- l3 %>% 
    filter(out == "prs") %>%
    group_by(h2_rare, what, h2, rho, gen) %>%
    summarise(r=median(r^2)) %>%
    mutate(
        rholab = paste("rho ==", rho), h2lab = paste("'Rare '* h^2 ==", h2_rare)
    )

temp %>%
    ggplot(., aes(x=as.factor(gen), y=r, colour=as.factor(what))) + geom_line(aes(group=what)) +
    facet_grid(h2lab ~ rholab, labeller = label_parsed) +
    labs(x="Generation", y=expression("Median" ~ R^2 ~ "between variant and PRS"), colour="Variant") +
    scale_colour_brewer(type="qual") +
    theme_minimal()
ggsave(file="results/run_sim_rare_causal.pdf", width=7, height=7)


load(file=here("results/run_sim4.RData"))
temp <- l4 %>% 
    filter(out == "yres") %>%
    group_by(h2_rare, what, h2, rho, gen, prs_known) %>%
    summarise(r=median(r^2)) %>%
    mutate(
        rholab = paste("rho ==", rho), h2lab = paste("'Rare '* h^2 ==", h2_rare)
    )

p1 <- temp %>%
    ggplot(., aes(x=gen, y=r)) + 
    geom_point(aes(colour=as.factor(prs_known))) +
    geom_smooth(aes(colour=as.factor(prs_known)), se=FALSE) +
    facet_grid(. ~ what) +
    labs(x="Generation", y=expression("Median" ~ R^2 ~ "between variant and trait adjusted for PRS"), linetype="Variant", colour="% h2 explained by PRS") +
    scale_colour_brewer(type="seq") +
    theme_minimal() +
    theme(axis.text.x= element_text(angle=90, hjust=0.5))
ggsave(p1, file="results/run_prs_adjustment.pdf", width=14)
