
# Print individual chapter
bookdown::preview_chapter("index.Rmd", "bookdown::pdf_book")
bookdown::preview_chapter("02-intro.Rmd", "bookdown::pdf_book")
bookdown::preview_chapter("03-data_exploration.Rmd", "bookdown::pdf_book")
bookdown::preview_chapter("04-bayes_general_model.Rmd", "bookdown::pdf_book")
bookdown::preview_chapter("05-bayes_poisson_glm.Rmd", "bookdown::pdf_book")
bookdown::preview_chapter("06-bayes_nb_glm.Rmd", "bookdown::pdf_book")
bookdown::preview_chapter("07-bayes_nb_bernoulli_glm.Rmd", "bookdown::pdf_book")
bookdown::preview_chapter("08-bayes-gamma-gml.Rmd", "bookdown::pdf_book")
bookdown::preview_chapter("09-bayes-implement-assess.Rmd", "bookdown::pdf_book")
bookdown::preview_chapter("10-bayes-coda.Rmd", "bookdown::pdf_book")

# Print book to PDF
bookdown::render_book("index.Rmd", "bookdown::pdf_book")

