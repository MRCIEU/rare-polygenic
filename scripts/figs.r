library(dplyr)
library(ggplot2)
library(here)

load(file=here("results/run_sim.RData"))

l %>% 
    filter(out == "prs") %>%
    group_by(gen, what, h2, rho) %>%
    summarise(r=median(r^2)) %>%
    mutate(gen=gen+1, what = case_when(what == "rare" ~ "Null rare variants", what == "common" ~ "Null common variants", what == "prs" ~ "Causal variants")) %>%
    mutate(rholab = paste("rho ==", rho), h2lab = paste("h^2 ==", h2)) %>%
    ggplot(., aes(x=as.factor(gen), y=r, colour=as.factor(what))) + geom_line(aes(group=what)) +
    facet_grid(h2 ~ rho, labeller = label_both) +
    labs(x="Generation", y="Median R^2 between variant and PRS", colour="Variant") +
    scale_colour_brewer(type="qual")
ggsave(here("results/run_sim_prs.pdf"), width=7, height=7)

l %>% 
    filter(out == "y") %>%
    group_by(gen, what, h2, rho) %>%
    summarise(r=median(r^2)) %>%
    mutate(gen=gen+1, what = case_when(what == "rare" ~ "Null rare variants", what == "common" ~ "Null common variants", what == "prs" ~ "Causal variants")) %>%
    mutate(rholab = paste("rho ==", rho), h2lab = paste("h^2 ==", h2)) %>%
    ggplot(., aes(x=as.factor(gen), y=r, colour=as.factor(what))) + geom_line(aes(group=what)) +
    facet_grid(h2lab ~ rholab, labeller = label_parsed) +
    labs(x="Generation", y=expression("Median" ~ R^2 ~ "between variant and trait"), colour="Variant") +
    scale_colour_brewer(type="qual")
ggsave(here("results/run_sim_trait.pdf"), width=7, height=7)

