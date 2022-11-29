all: build

build:
	(cd ./src; $(MAKE))

clean:
	(cd ./src; $(MAKE) clean)