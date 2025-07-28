# Название исполняемого файла
TARGET = lab3

# Все .cpp в текущей директории
SRCS := main.cpp func.cpp

# Цель по умолчанию
all: $(TARGET)

# Сборка: все .cpp компилируются сразу в один бинарник
$(TARGET): $(SRCS)
	g++ $(SRCS) -o $(TARGET)

# Построение графиков
graphs:
	python3 plt.py

report :
	cd doc && \
	rm -rf latex html && \
	doxygen Doxyfile && \
	cd latex && \
	make

# Запуск
run: all
	./$(TARGET)

# Полный запуск
full: clean-all all run graphs report

# Очистка результатов
clean-result:
	rm -f times.csv

# Очистка исполняемого файла
clean:
	rm -f $(TARGET)

# Полная очистка
clean-all: clean-result clean
