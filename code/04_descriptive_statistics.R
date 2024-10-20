# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hpv/code/00_package.R",echo = T)

library(ggrepel)
# Read aggregated data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("hpv/clean_data/aggregated_data.rda")

transorm_fun <- function(x){
  x %>% 
    mutate(
    name = case_when(
      name == "hpv16" ~ "Positive HPV-16",
      name == "hpv18" ~ "Positive HPV-18",
      name == "hpvneg" ~ "Negative to all genotypes",
      name == "hpvother" ~ "Positive HPV-other"
    ),
    name = factor(name, levels = c("Positive HPV-16",
                                   "Positive HPV-18",
                                   "Positive HPV-other",
                                   "Negative to all genotypes")),
    type = case_when(type == "ind" ~ "Opportunistic", 
                     type == "orga" ~ "Organised",
                     TRUE ~ "All")
  )
}

formating_text <- function(x,
                           digits_choice = 0) {
  format(
    round(x,
          digits = digits_choice
    ),
    small.mark = ".",
    big.mark = ",",
    nsmall = digits_choice
  )
}

size_font <- 20

data_bru <- data_bru %>% st_drop_geometry()
# N global ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sample_total <- data_bru %>%
  rbind(data_bru %>% 
          mutate(type = "All")) %>% 
  group_by(type,year) %>%
  summarise(N  = sum(hpvneg) + sum(hpvpos))  %>% 
  pivot_longer(cols = c(-type,-year),
               names_to = "name",
               values_to = "value") %>% 
  mutate(
    type = case_when(type == "ind" ~ "Opportunistic", 
                     type == "orga" ~ "Organised",
                     TRUE ~ "All")) %>% 
  ungroup() %>% 
  mutate(label = ifelse(type != "All", formating_text(value,0), NA)) %>% 
  ggplot(aes(x=year, y=value, color = type)) + 
  geom_point(position="dodge", stat="identity", size = 3)+ 
  geom_line(position="dodge", stat="identity")+
  labs(x = "Year",
       y = "Number of tests",
       color = "Type of screening: ") + 
  scale_color_manual(values = c("black", "red", "blue")) +
  scale_y_continuous(labels = scales::comma, n.breaks = 10) +
  theme_minimal(base_size = size_font)+
  theme(legend.position = "bottom") +
  geom_text(
    aes(label = label, 
        y = value -12500,
        color = type
    ),
    bg.color = "white",
    size = 5
  )


number_test_age <- data_bru %>%
  rbind(data_bru %>% 
          mutate(type = "All")) %>% 
  group_by(type,age_class) %>%
  summarise(N  = sum(hpvneg) + sum(hpvpos))  %>% 
  pivot_longer(cols = c(-type,-age_class),
               names_to = "name",
               values_to = "value") %>% 
  mutate(
    type = case_when(type == "ind" ~ "Opportunistic", 
                     type == "orga" ~ "Organised",
                     TRUE ~ "All")) %>% 
  ungroup() %>% 
  mutate(label = ifelse(type != "All", formating_text(value,0), NA)) %>% 
  filter(type != "All") %>% 
  ggplot(aes(x=age_class, y=value, fill = type, color = type)) + 
  geom_bar(position="dodge", stat="identity", alpha = 0.5)+
  labs(x = "Age",
       y = "Number of tests",
       color = "Type of screening: ",
       fill = "Type of screening: ") + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", 'blue')) +
  scale_y_continuous(labels = scales::comma, n.breaks = 10) +
  theme_minimal(base_size = size_font)+
  theme(legend.position = "bottom")  +
  theme(axis.text.x = element_text(angle = +90, vjust = 0.5, hjust=1)) +
  facet_wrap(~type, scale = "free")+
  theme(legend.position = "none") 


# % year-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
descriptive_plot_year <- data_bru %>%
  rbind(data_bru %>% 
          mutate(type = "All")) 

descriptive_plot_year <- descriptive_plot_year %>% 
  rbind(descriptive_plot_year %>% 
          mutate(year = "All")) %>% 
  mutate(year = factor(year, levels = c("All", 
                                        "2020",
                                        "2021",
                                        "2022",
                                        "2023")))

descriptive_plot_year <- descriptive_plot_year %>% 
  group_by(type,year) %>%
  summarise(hpv16 = sum(hpv16),
            hpv18 = sum(hpv18),
            hpvother = sum(hpvother),
            hpvneg  = sum(hpvneg),
            hpvpos = sum(hpvpos)) %>% 
  mutate(N = hpvneg + hpvpos) %>% 
  mutate_at(
    vars("hpv16",
         "hpv18",
         "hpvother",
         "hpvneg",
         "hpvpos"),
    ~./N) %>% 
  pivot_longer(cols = c(-type,-year),
               names_to = "name",
               values_to = "value") %>% 
  transorm_fun(.) %>% 
  filter(!is.na(name)) %>% 
  mutate(year = factor(year)) %>% 
  ggplot(aes(x=type, y=value, fill = year)) + 
  geom_bar(position="dodge", stat="identity", alpha = 0.5)+
  labs(x = "Type of screening",
       y = "Proportion of tests",
       fill = "Year: ") + 
  scale_fill_manual(values = c("black", "pink", "orange", "cyan", "purple")) +
  scale_y_continuous(n.breaks = 5,label = scales::percent) +
  theme_minimal(base_size = size_font)+
  theme(legend.position = "bottom") +
  facet_wrap(~name, scale = "free")

# Descriptive plot age ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
descriptive_plot_age <- data_bru %>%
  group_by(age_class,type) %>%
  summarise(hpv16 = sum(hpv16),
            hpv18 = sum(hpv18),
            hpvother = sum(hpvother),
            hpvneg  = sum(hpvneg),
            hpvpos = sum(hpvpos)) %>% 
  mutate(N = hpvneg + hpvpos) %>% 
  mutate_at(
    vars("hpv16",
         "hpv18",
         "hpvother",
         "hpvneg",
         "hpvpos"),
    ~./N) %>% 
  pivot_longer(cols = c(-age_class, -type),
               names_to = "name",
               values_to = "value") %>% 
  transorm_fun(.) %>% 
  filter(!is.na(name)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(aes(y = value, x = age_class, color = type), size = 3) +
  scale_color_manual(values = c("red", "blue")) +
  facet_wrap( ~ name, scale = "free") +
  scale_y_continuous(label = scales::percent) +
  labs(x = "Age class",
       y = "Proportion of tests",
       color = 'Type of screening:') +
  theme_minimal(base_size = size_font)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.position = "bottom") 

# N time ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tmp <- data_bru %>% 
  rbind(data_bru %>% mutate(type = "All")) %>% 
  group_by(type,week, idtime) %>% 
  summarise(N = sum(N)) %>%
  ungroup() %>% 
  dplyr::select(-idtime) %>% 
  right_join(dfdate, by = "week") %>% 
  filter(!is.na(type)) %>% 
  mutate(type = case_when(type == "All" ~ "All",
                          type == "ind" ~ "Opportunistic",
                          type == "orga" ~ "Organised"))

label <- tmp %>% 
  filter(row_number() %in% seq(1,nrow(tmp),50) |
          week %in% c("2021-01", "2022-01", "2023-01")) %>% 
  filter(type == "All")

plot_ntest <- tmp %>% 
  ggplot() +
  geom_line(aes(x = idtime,
                y = N), linewidth = 1 )+
  labs(x = "Time (Year-Week)",
       y = "Number of tests") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = label$idtime,
                     label = label$week)  +
  theme_minimal(base_size = size_font) +
  facet_wrap(~type, scale = "free_y", ncol = 1) +
  theme(axis.text.x = element_text(angle = +90, vjust = 0.5, hjust=1)) 
  


# Positivity --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------<
tmp <- data_bru %>% 
  rbind(data_bru %>% mutate(type = "All")) %>% 
  group_by(type,week, idtime) %>%
  summarise(hpv16 = sum(hpv16),
            hpv18 = sum(hpv18),
            hpvother = sum(hpvother),
            hpvneg  = sum(hpvneg),
            hpvpos = sum(hpvpos)) %>% 
  mutate(N = hpvneg + hpvpos) %>% 
  mutate_at(
    vars("hpv16",
         "hpv18",
         "hpvother",
         "hpvneg",
         "hpvpos"),
    ~./N) %>% 
  pivot_longer(cols = c(-type, -week, -idtime),
               names_to = "name",
               values_to = "value") %>% 
  transorm_fun(.) %>% 
  filter(!is.na(name)) %>% 
  ungroup() 


label <- tmp %>% 
  filter(type == "All") %>% 
  filter(name == "Positive HPV-16") %>% 
  filter(row_number() %in% seq(1,nrow(tmp),50) |
           week %in% c("2021-01", "2022-01", "2023-01"))


plot_pcttest <- tmp %>% 
  ggplot() +
  geom_point(aes(x = idtime,
                y = value,
                color = type))+
  labs(x = "Time (Year-Week)",
       y = "Number of tests",
       color = 'Type of tests: ') +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = label$idtime,
                     label = label$week)  +
  theme_minimal(base_size = size_font) +
  facet_wrap(~name, scale = "free_y") +
  theme(axis.text.x = element_text(angle = +90, vjust = 0.5, hjust=1)) 


# Save everything ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
out <- list(
  "descriptive_plot_sample" = sample_total,
  "descriptive_plot_year" = descriptive_plot_year,
  "number_test_age" = number_test_age,
  "descriptive_plot_age" = descriptive_plot_age,
  "plot_ntest" = plot_ntest,
  "plot_pcttest" = plot_pcttest,
  "transorm_fun" = transorm_fun
)

saveRDS(out,
        "hpv/clean_data/descriptive_figures.RDS")